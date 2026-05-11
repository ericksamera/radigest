#!/usr/bin/env python3
"""Run radigest across many enzyme pairs and write one JSON per pair.

This script keeps radigest as the digest/scoring engine. It is intended for
large ddRAD/GBS design screens where only JSON summaries are needed initially.
"""

from __future__ import annotations

import argparse
import itertools
import re
import shlex
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class PairJob:
    enzyme_a: str
    enzyme_b: str
    tag: str
    json_path: Path
    log_path: Path
    command: list[str]


@dataclass(frozen=True)
class PairResult:
    tag: str
    returncode: int
    json_path: Path
    log_path: Path


SAFE_TAG = re.compile(r"[^A-Za-z0-9_.-]+")


def read_enzyme_names(value: str) -> list[str]:
    path = Path(value)
    if path.exists():
        names: list[str] = []
        for line in path.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            names.extend(part.strip() for part in line.split(",") if part.strip())
    else:
        names = [part.strip() for part in value.split(",") if part.strip()]

    deduped: list[str] = []
    seen: set[str] = set()
    for name in names:
        if name not in seen:
            deduped.append(name)
            seen.add(name)
    if len(deduped) < 2:
        raise SystemExit("Need at least two candidate enzymes.")
    return deduped


def pair_tag(a: str, b: str) -> str:
    raw = f"{a}__{b}"
    return SAFE_TAG.sub("_", raw)


def build_jobs(args: argparse.Namespace, enzymes: list[str]) -> list[PairJob]:
    radigest_cmd = shlex.split(args.radigest)
    json_dir = args.out_dir / "json"
    log_dir = args.out_dir / "logs"
    json_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    pairs = itertools.combinations(enzymes, 2)
    jobs: list[PairJob] = []
    for enzyme_a, enzyme_b in pairs:
        tag = pair_tag(enzyme_a, enzyme_b)
        json_path = json_dir / f"{tag}.json"
        log_path = log_dir / f"{tag}.log"
        if json_path.exists() and not args.force:
            continue

        command = [
            *radigest_cmd,
            "-fasta",
            str(args.fasta),
            "-enzymes",
            f"{enzyme_a},{enzyme_b}",
            "-min",
            str(args.min),
            "-max",
            str(args.max),
            "-score-min",
            str(args.score_min),
            "-score-max",
            str(args.score_max),
            "-size-model",
            args.size_model,
            "-size-mean",
            str(args.size_mean),
            "-size-sd",
            str(args.size_sd),
            "-size-edge-sd",
            str(args.size_edge_sd),
            "-threads",
            str(args.radigest_threads),
            "-fragments-tsv",
            "",
            "-gff",
            args.gff_sink,
            "-json",
            str(json_path),
        ]
        if args.allow_same:
            command.append("-allow-same")
        for extra in args.extra_arg:
            command.extend(shlex.split(extra))

        jobs.append(PairJob(enzyme_a, enzyme_b, tag, json_path, log_path, command))
        if args.max_pairs is not None and len(jobs) >= args.max_pairs:
            break
    return jobs


def run_job(job: PairJob, dry_run: bool) -> PairResult:
    if dry_run:
        print(shlex.join(job.command))
        return PairResult(job.tag, 0, job.json_path, job.log_path)

    proc = subprocess.run(job.command, check=False, capture_output=True, text=True)
    with job.log_path.open("w") as log:
        log.write("$ " + shlex.join(job.command) + "\n\n")
        log.write("# stdout\n")
        log.write(proc.stdout)
        if proc.stdout and not proc.stdout.endswith("\n"):
            log.write("\n")
        log.write("\n# stderr\n")
        log.write(proc.stderr)
        if proc.stderr and not proc.stderr.endswith("\n"):
            log.write("\n")
    return PairResult(job.tag, proc.returncode, job.json_path, job.log_path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run radigest for all unique pairs in a candidate enzyme set."
    )
    parser.add_argument("--fasta", type=Path, required=True, help="reference FASTA")
    parser.add_argument(
        "--enzymes",
        required=True,
        help="comma-separated enzymes or a file with one/comma-separated enzymes per line",
    )
    parser.add_argument(
        "--radigest",
        default="go run ./cmd/radigest",
        help="radigest command, e.g. 'radigest' or 'go run ./cmd/radigest'",
    )
    parser.add_argument("--out-dir", type=Path, default=Path("pair_screen"))
    parser.add_argument("--min", type=int, default=300)
    parser.add_argument("--max", type=int, default=600)
    parser.add_argument("--score-min", type=int, default=1)
    parser.add_argument("--score-max", type=int, default=2000)
    parser.add_argument("--size-model", default="normal")
    parser.add_argument("--size-mean", type=float, default=275.0)
    parser.add_argument("--size-sd", type=float, default=85.0)
    parser.add_argument("--size-edge-sd", type=float, default=25.0)
    parser.add_argument("--radigest-threads", type=int, default=1)
    parser.add_argument("--jobs", type=int, default=1, help="parallel radigest processes")
    parser.add_argument("--gff-sink", default="/dev/null")
    parser.add_argument("--allow-same", action="store_true")
    parser.add_argument("--force", action="store_true", help="overwrite existing JSONs")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--max-pairs", type=int, default=None)
    parser.add_argument(
        "--extra-arg",
        action="append",
        default=[],
        help="extra radigest argument string; may be repeated",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    enzymes = read_enzyme_names(args.enzymes)
    jobs = build_jobs(args, enzymes)
    print(f"candidate_enzymes\t{len(enzymes)}", file=sys.stderr)
    print(f"pair_jobs\t{len(jobs)}", file=sys.stderr)
    if not jobs:
        return 0

    failures: list[PairResult] = []
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        futures = [pool.submit(run_job, job, args.dry_run) for job in jobs]
        for future in as_completed(futures):
            result = future.result()
            status = "ok" if result.returncode == 0 else f"exit={result.returncode}"
            print(f"{result.tag}\t{status}\t{result.json_path}\t{result.log_path}", file=sys.stderr)
            if result.returncode != 0:
                failures.append(result)

    if failures:
        print("failed_pairs", len(failures), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
