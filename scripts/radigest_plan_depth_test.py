#!/usr/bin/env python3
from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).with_name("radigest-plan-depth")
HEADER = [
    "rank",
    "enzyme_a",
    "enzyme_b",
    "score",
    "weighted_genome_pct",
    "raw_genome_pct",
    "target_delta_pct",
    "abs_target_delta_pct",
    "weighted_bases",
    "weighted_fragments",
    "raw_bases_in_window",
    "raw_fragments_in_window",
    "mean_weighted_length",
    "json_path",
]


def run_plan(
    row: dict[str, str],
    *extra: str,
    budget_args: tuple[str, str] = ("--lane-read-pairs", "300M"),
) -> dict[str, str]:
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "ranked.tsv"
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=HEADER, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerow(row)
        cmd = [
            sys.executable,
            str(SCRIPT),
            str(path),
            "--read-layout",
            "pe",
            "--read-length",
            "150",
            *budget_args,
            "--usable-read-fraction",
            "0.8",
            "--desired-depth",
            "10",
            *extra,
        ]
        proc = subprocess.run(cmd, check=True, text=True, capture_output=True)
        rows = list(csv.DictReader(proc.stdout.splitlines(), delimiter="\t"))
        assert len(rows) == 1
        return rows[0]


class PlanDepthTest(unittest.TestCase):
    def base_row(self) -> dict[str, str]:
        return {
            "rank": "1",
            "enzyme_a": "EcoRI",
            "enzyme_b": "MseI",
            "score": "2.5",
            "weighted_genome_pct": "2.5",
            "raw_genome_pct": "3.0",
            "target_delta_pct": "",
            "abs_target_delta_pct": "",
            "weighted_bases": "25000000",
            "weighted_fragments": "100000",
            "raw_bases_in_window": "30000000",
            "raw_fragments_in_window": "120000",
            "mean_weighted_length": "250",
            "json_path": "pair.json",
        }

    def test_samples_depth_and_target_are_solved(self) -> None:
        out = run_plan(self.base_row(), "--samples", "96", "--target-genome-pct", "1.25")
        self.assertEqual(out["read_pairs_per_sample"], "2500000.000000")
        self.assertEqual(out["mean_depth_full_target"], "25.000000")
        self.assertEqual(out["full_target_passes_depth"], "true")
        self.assertEqual(out["budget_supported_genome_pct"], "2.500000")
        self.assertEqual(out["required_pairs_per_sample_full_target"], "1000000.000000")
        self.assertEqual(out["max_samples_per_lane_full_target"], "240")
        self.assertEqual(out["target_possible"], "true")
        self.assertEqual(out["target_fraction_of_weighted_target"], "0.500000")
        self.assertEqual(out["target_weighted_fragments"], "50000.000000")
        self.assertEqual(out["max_samples_per_lane_for_target"], "480")
        self.assertEqual(out["mean_insert_category"], "mean_lt_2_read_lengths_overlap_risk")

    def test_capacity_mode_without_samples_reports_max_samples(self) -> None:
        out = run_plan(self.base_row(), "--target-genome-pct", "1.25")
        self.assertEqual(out["samples"], "")
        self.assertEqual(out["read_pairs_per_sample"], "")
        self.assertEqual(out["max_samples_total_for_target"], "480")
        self.assertEqual(out["target_possible"], "true")

    def test_target_impossible_leaves_target_capacity_blank(self) -> None:
        out = run_plan(self.base_row(), "--target-genome-pct", "3.0")
        self.assertEqual(out["target_possible"], "false")
        self.assertEqual(out["target_weighted_fragments"], "")
        self.assertEqual(out["required_pairs_per_sample_for_target"], "")
        self.assertEqual(out["max_samples_total_for_target"], "")

    def test_flowcell_read_pairs_alias_matches_single_lane_budget(self) -> None:
        out = run_plan(
            self.base_row(),
            "--samples",
            "96",
            "--target-genome-pct",
            "1.25",
            budget_args=("--flowcell-read-pairs", "300M"),
        )
        self.assertEqual(out["lane_read_pairs"], "300000000.000000")
        self.assertEqual(out["lanes"], "1")
        self.assertEqual(out["max_samples_total_for_target"], "480")

    def test_budget_supported_genome_pct_scales_when_depth_is_short(self) -> None:
        out = run_plan(self.base_row(), "--samples", "600")
        self.assertEqual(out["mean_depth_full_target"], "4.000000")
        self.assertEqual(out["full_target_passes_depth"], "false")
        self.assertEqual(out["budget_supported_genome_pct"], "1.000000")
        self.assertEqual(out["budget_supported_weighted_bases"], "10000000.000000")

    def test_genome_denominator_backfills_missing_weighted_genome_pct(self) -> None:
        row = self.base_row()
        row["weighted_genome_pct"] = ""
        out = run_plan(row, "--genome-bases", "1G", "--samples", "96")
        self.assertEqual(out["weighted_genome_pct"], "2.500000")
        self.assertEqual(out["budget_supported_genome_pct"], "2.500000")


if __name__ == "__main__":
    unittest.main()
