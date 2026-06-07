# radigest

[![CI](https://img.shields.io/github/actions/workflow/status/ericksamera/radigest/ci.yml?branch=main&label=ci)](https://github.com/ericksamera/radigest/actions/workflows/ci.yml)
[![DOI](http://zenodo.org/badge/979818941.svg)](https://doi.org/10.5281/zenodo.20176743)
[![Go](https://img.shields.io/badge/go-%3E=%201.22-blue)](https://golang.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)

Fast in-silico restriction digest for genomics. Give it a reference FASTA (plain or `.gz`) or synthesize one on the fly; pass one or two enzymes; it scans, size-selects, and reports a JSON run summary by default. Optional GFF3, BED, TSV, and FASTA fragment exports are available for GBS/ddRAD, probe design, visualization, and downstream modeling. Output order is deterministic even with multithreading.

---

## Features

- **Single or double digest.** Double-digest keeps **adjacent AB/BA** by default; enable **AA/BB** too with `-allow-same`. Single-digest uses consecutive A cuts. Terminal chromosome/contig-end fragments are omitted by default; keep them with `-include-ends`.
- **IUPAC & cut offsets.** Sites accept degenerate codes; the cut index comes from `^` in the site (or mid-site if missing). `-strict-cuts` makes missing carets an error.
- **Robust FASTA I/O.** Read from a path or `-` (STDIN), auto-detect `.gz`, normalize case, and **trim CRLF**. `N` in the **reference** does **not** match any site.
- **Synthetic genomes.** Generate a single-chromosome genome named `chr1` with `-sim-len`, `-sim-gc`, `-sim-seed` and digest it directly—no FASTA on disk needed.
- **Clean outputs.** With no output flags, radigest writes a JSON summary to STDOUT. GFF3, BED, fragment FASTA, and per-fragment TSV are opt-in artifact outputs. GFF3 uses `ID=<chr>_<n>;Length=<bp>`; TSV records insert length, hard-keep status, and size-selection weight. Coordinates are **1-based closed** in GFF and **0-based half-open** in BED, TSV, and FASTA metadata.
- **Size-selection scoring.** Keep the hard `-min/-max` window for hard-kept outputs while assigning per-fragment recovery weights with `hard`, `normal`, `triangular`, or `soft-window` models over an optional broader score range.
- **Streaming fragment export.** The CLI streams digest fragments to the collector instead of materializing every kept fragment for a chromosome before writing.

---

## Build and install

```bash
# Public tool surface: radigest, radigest-design, radigest-fit-size-model
make build
make install

# Developer/helper surface: cached screening, benchmarks, and legacy helper scripts
make build-dev
make install-dev
```

The default install target keeps user-facing commands small and intentional. `make install-dev` additionally installs `radigest-screen-pairs-cached`, `radigest-bench-screen-cached`, `radigest-screen-pairs`, `radigest-rank-pairs`, and `radigest-plan-depth` for development and diagnostic workflows.

---

## Quick start

```bash
# Default: write a run summary JSON document to stdout
radigest -fasta ref.fa -enzymes EcoRI,MseI

# Single digest (EcoRI) → GFF and BED files
radigest -fasta ref.fa -enzymes EcoRI -gff fragments.gff3 -bed fragments.bed

# Double digest with size selection + JSON summary file
radigest -fasta ref.fa -enzymes EcoRI,MseI -min 100 -max 800 -json run.json

# Also write FASTA sequences for the hard-kept fragments
radigest -fasta ref.fa -enzymes EcoRI,MseI -min 100 -max 800 \
  -gff fragments.gff3 -bed fragments.bed -fragments-fasta fragments.fa

# ddRAD-style soft-window scoring with broad per-fragment TSV for downstream modeling
radigest -fasta ref.fa -enzymes PstI,MspI -min 250 -max 500 -score-min 1 -score-max 1000 \
  -size-model soft-window -size-edge-sd 25 -fragments-tsv fragments.tsv -json run.json

# Double digest but ALSO keep AA/BB neighbors
radigest -fasta ref.fa -enzymes EcoRI,MseI -allow-same -gff fragments.gff3

# Include terminal fragments from chromosome/contig ends to the nearest cut
radigest -fasta ref.fa -enzymes EcoRI,MseI -include-ends -gff fragments.gff3

# Simulate a 10 Mb genome at 42% GC and digest, reporting JSON to stdout
radigest -sim-len 10000000 -sim-gc 0.42 -sim-seed 123 -enzymes EcoRI,MseI
```

---

## CLI (most used)

- `-fasta <path|->` — reference FASTA; `-` = STDIN; `.gz` auto-detected.
- `-enzymes E1[,E2]` — one (A) or two (A,B) only. In double-digest, AB/BA by default.
- `-min/-max` — hard-selected insert-length window used for hard-kept artifact outputs and `hard_kept` in TSV (**default min=1**).
- `-score-min/-score-max` — broader insert-length range used for weighted size-selection stats and emitted to `-fragments-tsv` when TSV output is enabled; defaults to `-min/-max`.
- `-size-model hard|normal|triangular|soft-window` — per-fragment size-selection weight model (**default `hard`**).
- `-size-mean`, `-size-sd`, `-size-edge-sd` — parameters for `normal`, `triangular`, and `soft-window` scoring.
- `-json <path|->` — write a run summary JSON document with counts, bases, per-chromosome stats, and size-selection weighted stats; `-` = STDOUT. When no output flags are set, radigest writes this JSON summary to STDOUT by default.
- `-gff <path|->` — optional GFF3 output for hard-kept fragments; `-` = STDOUT; empty string disables (**default disabled**).
- `-bed <path|->` — optional BED6 output for hard-kept fragments; `-` = STDOUT; empty string disables (**default disabled**).
- `-fragments-fasta <path|->` — optional FASTA sequences for hard-kept fragments, using the same saved fragment set and ordinals as GFF; `-` = STDOUT; empty string disables (**default disabled**).
- `-fragments-tsv <path|->` — optional per-fragment TSV for score-range fragments; `-` = STDOUT; empty string disables (**default disabled**).
- `-threads <n>` — positive worker count; `-v`, `-version`, `-list-enzymes`.
- **Simulation:** `-sim-len <bp>`, `-sim-gc <0..1>` (invalid values error), `-sim-seed <int>` (emits a single `chr1`).
- **Modes:** `-allow-same` (keep AA/BB in double-digest), `-include-ends` (also emit terminal chromosome/contig-end fragments), `-strict-cuts` (error if a site lacks `^` and would otherwise fall back to mid-site).

---

## Scope and limitations

radigest is a deterministic sequence-level model. It identifies recognition sites and cut coordinates from the reference sequence only. It does **not** model methylation sensitivity, partial digestion, star activity, enzyme efficiency, buffer compatibility, or empirical digestion rates. Enzymes with the same recognition motif and cut coordinate are treated identically by the digest logic even when their wet-lab behavior can differ under methylation or assay conditions.

---

## Outputs

If no output flags are set, radigest writes a run summary JSON document to STDOUT. If any output flag is set, radigest writes exactly the requested outputs. GFF3, BED, TSV, and FASTA are disabled unless explicitly requested. At most one active output may target STDOUT.

### JSON summary

```json
{
  "schema_version": 1,
  "radigest_version": "v0.1.0",
  "command": [
    "radigest",
    "-fasta",
    "ref.fa",
    "-enzymes",
    "EcoRI,MseI",
    "-min",
    "100",
    "-max",
    "800",
    "-score-min",
    "1",
    "-score-max",
    "1000",
    "-size-model",
    "soft-window",
    "-size-edge-sd",
    "25",
    "-json",
    "run.json",
    "-gff",
    "fragments.gff3",
    "-bed",
    "fragments.bed",
    "-fragments-tsv",
    "fragments.tsv",
    "-fragments-fasta",
    "fragments.fa",
    "-threads",
    "8"
  ],
  "input": {
    "source": "fasta",
    "fasta": "ref.fa"
  },
  "parameters": {
    "min_length": 100,
    "max_length": 800,
    "score_min": 1,
    "score_max": 1000,
    "size_model": "soft-window",
    "size_edge_sd": 25,
    "threads": 8,
    "allow_same": false,
    "strict_cuts": false,
    "include_ends": false
  },
  "outputs": {
    "json": "run.json",
    "gff": "fragments.gff3",
    "bed": "fragments.bed",
    "fragments_tsv": "fragments.tsv",
    "fragments_fasta": "fragments.fa"
  },
  "warnings": [],
  "enzymes": ["EcoRI", "MseI"],
  "min_length": 100,
  "max_length": 800,
  "gff": "fragments.gff3",
  "bed": "fragments.bed",
  "fragments_tsv": "fragments.tsv",
  "fragments_fasta": "fragments.fa",
  "size_selection": {
    "model": "soft-window",
    "score_min": 1,
    "score_max": 1000,
    "edge_sd": 25,
    "raw_fragments_scored": 234567,
    "raw_bases_scored": 91234567,
    "raw_fragments_in_window": 123456,
    "raw_bases_in_window": 42100000,
    "weighted_fragments": 98234.7,
    "weighted_bases": 33100000.5,
    "mean_weighted_length": 336.9
  },
  "total_fragments": 123456,
  "total_bases": 7891011,
  "per_chromosome": { "chr1": { "fragments": 23456, "bases": 3456789 } }
}
```

`schema_version`, `radigest_version`, `command`, `input`, `parameters`, `outputs`, and `warnings` provide a stable provenance header for downstream tools. For simulated input, the `input` block records both `sim_seed_requested` and `sim_seed_resolved`; the resolved seed reproduces runs where `-sim-seed 0` requested a time-based seed. The top-level `gff`, `bed`, `fragments_tsv`, and `fragments_fasta` fields are retained for compatibility and are present only when those artifact outputs are enabled.

### GFF3

```
##gff-version 3
chr1	radigest	fragment	<start>	<end>	.	+	.	ID=chr1_1;Length=123
```

`start/end` are **1-based closed**; `Length` is `end - start + 1`. Ordering is deterministic per chromosome.

### BED

When `-bed` is set, radigest writes hard-kept fragments as BED6 records with no header:

```text
chr1	10422	10731	chr1_1	0	+
```

`chromStart/chromEnd` are **0-based half-open**. The fourth column is the same per-chromosome fragment ordinal used by GFF and fragment FASTA; the score is `0` and strand is `+`.

### Fragment FASTA

When `-fragments-fasta` is set, radigest writes FASTA records for the hard-kept fragments:

```text
>chr1_1 chrom=chr1 start0=10422 end0=10731 length=309
AATT...
```

The fragment ID uses the same per-chromosome ordinal as GFF. Header coordinates are 0-based half-open. Use `-min 0 -max <large>` to emit every internal digest fragment that radigest would otherwise keep under the selected digest mode; terminal contig-end fragments still require `-include-ends`.

### Fragment TSV

When `-fragments-tsv` is set, radigest writes a per-fragment TSV for fragments in the score range:

```
chrom	start0	end0	length	hard_kept	size_weight
chr1	10422	10731	309	true	0.982143
chr1	18831	18922	91	false	0.014221
```

`hard_kept` is true when the insert length is inside `-min/-max`. `size_weight` is the selected size model evaluated on insert length only.

---

## Pair-screening helper scripts

The Go CLI intentionally stays focused on digesting and fragment scoring. For larger ddRAD/GBS design screens, use the helper scripts in `scripts/` to run many enzyme pairs and rank the resulting JSON summaries.

Create a candidate enzyme list:

```bash
cat > candidate_enzymes.txt <<'EOF2'
EcoRI
MseI
PstI
MspI
ApeKI
NlaIII
MluCI
BfaI
EOF2
```

Run every unique pair using an empirically calibrated size model. This example uses the sockeye ddRAD profile fitted from observed TLENs, `normal(mean=275, sd=85)`:

```bash
scripts/radigest-screen-pairs \
  --fasta ref.fa \
  --enzymes candidate_enzymes.txt \
  --min 300 \
  --max 600 \
  --score-min 1 \
  --score-max 2000 \
  --size-model normal \
  --size-mean 275 \
  --size-sd 85 \
  --jobs 2 \
  --radigest-threads 2 \
  --out-dir pair_screen
```

The screen writes one JSON summary per pair under `pair_screen/json/` and logs under `pair_screen/logs/`. It requests only JSON summaries, so GFF, BED, TSV, and FASTA artifact outputs are omitted during initial screening.

For diagnostic timing of the cached pair-screen engine itself, use the compiled benchmark binary. It reports separate cut-index build, pair-scoring, JSON-marshalling, and JSON-writing phases as TSV:

```bash
radigest-bench-screen-cached \
  --fasta ref.fa \
  --enzymes candidate_enzymes.txt \
  --min 300 \
  --max 600 \
  --score-min 1 \
  --score-max 2000 \
  --size-model hard \
  --jobs 4 \
  --build-workers 4 \
  --runs 3 \
  --reuse-index \
  --output-mode none
```

Use `--build-workers` to control parallel cut-index construction; when omitted, cached screening defaults it from `--jobs` and then `--threads`. Set `--build-workers 1` to reproduce the older serial cache-build behavior. Use `--output-mode marshal` to add JSON encoding cost without filesystem writes, or `--output-mode write --out-dir pair_screen_bench --force` to include per-pair JSON file writes.

Rank pairs by weighted bases, or by genome percentage if a FASTA denominator is provided:

```bash
# Rank by weighted recovered insert-bases
scripts/radigest-rank-pairs 'pair_screen/json/*.json' \
  --objective weighted-bases \
  --out pair_screen/ranked_pairs.tsv

# Rank by weighted genome percentage using non-N reference bases as denominator
scripts/radigest-rank-pairs 'pair_screen/json/*.json' \
  --fasta ref.fa \
  --objective weighted-genome-pct \
  --out pair_screen/ranked_pairs.genome_pct.tsv

# Find pairs closest to a target weighted genome percentage
scripts/radigest-rank-pairs 'pair_screen/json/*.json' \
  --fasta ref.fa \
  --objective closest-target \
  --target-genome-pct 2.5 \
  --out pair_screen/ranked_pairs.closest_2.5pct.tsv

# Reuse a known non-N genome denominator without rereading the FASTA
GENOME_BASES=2643888753
scripts/radigest-rank-pairs 'pair_screen/json/*.json' \
  --genome-bases "$GENOME_BASES" \
  --objective closest-target \
  --target-genome-pct 1.5 \
  --out pair_screen/ranked_pairs.closest_1.5pct.tsv
```

The ranked TSV includes `weighted_bases`, `weighted_fragments`, `raw_bases_in_window`, `raw_fragments_in_window`, `mean_weighted_length`, and genome-percentage columns when a denominator is supplied.

Plan sequencing depth and multiplexing from the ranked TSV:

```bash
scripts/radigest-plan-depth pair_screen/ranked_pairs.genome_pct.tsv \
  --read-layout pe \
  --read-length 150 \
  --lane-read-pairs 300M \
  --lanes 1 \
  --usable-read-fraction 0.85 \
  --samples 96 \
  --desired-depth 10 \
  --target-genome-pct 2.5 \
  --out pair_screen/depth_plan.tsv

# MiSeq-style one-flowcell budget
scripts/radigest-plan-depth pair_screen/ranked_pairs.genome_pct.tsv \
  --read-layout pe \
  --read-length 150 \
  --flowcell-read-pairs 50M \
  --usable-read-fraction 0.85 \
  --samples 37 \
  --desired-depth 30 \
  --target-genome-pct 1.5 \
  --out pair_screen/depth_plan.flowcell.tsv
```

This first-pass planner treats `weighted_fragments` as the number of recovered loci competing for reads. `--desired-depth` is a mean read-pair depth per recovered locus, not a callable-depth distribution model or basewise WGS depth. It reports expected mean depth over the full weighted target, the genome percentage supported at the requested depth for a planned sample count, and the number of samples supportable per lane, flowcell, or run for a target genome percentage. It does not model per-locus depth dispersion.

---

## Inverse enzyme-pair design

Use `radigest-design` when the experimental design target is known and the enzyme pair is the unknown. The command screens candidate enzyme pairs with the cached cut-index engine, computes weighted recovered genome percentage, combines each pair with a sequencing budget, and ranks pairs by modeled fit to the requested genome fraction and mean read-pair depth.

```bash
radigest-design \
  --ref ref.fa \
  --pct 2.5 \
  --enzymes EcoRI,MseI,PstI,ApeKI,NlaIII,MspI \
  --samples 96 \
  --flowcell-read-pairs 300M \
  --depth 10 \
  --read-length 150
```

The concise aliases map to the explicit flags: `--ref` = `--fasta`, `--pct` = `--target-genome-pct`, and `--depth` = `--desired-depth`. The longer form remains available for scripts and reports:

```bash
radigest-design \
  --fasta ref.fa \
  --enzymes candidate_enzymes.txt \
  --target-genome-pct 2.5 \
  --coverage-tolerance-pct 0.25 \
  --desired-depth 10 \
  --samples 96 \
  --read-layout pe \
  --read-length 150 \
  --lane-read-pairs 300M \
  --lanes 1 \
  --usable-read-fraction 0.85 \
  --min 300 \
  --max 600 \
  --score-min 1 \
  --score-max 2000 \
  --size-model normal \
  --size-mean 275 \
  --size-sd 85 \
  --out-dir radigest_design
```

The command writes `design.summary.tsv`, `design.tsv`, and `design.json`. `design.summary.tsv` is the compact human-review table with rank, enzyme pair, feasibility, target and predicted genome percentage, target and predicted mean locus depth, read-pair budget, max samples, weighted fragments, mean insert, insert status, and `fit_score`. The full `design.tsv` keeps the machine-readable table with `feasible`, `decision_reason`, `fit_score`, `fit_loss`, `predicted_weighted_genome_pct`, `target_mean_locus_depth`, `predicted_mean_locus_depth`, depth shortfall, sequencing-budget columns, insert-size diagnostics, and cached-screening provenance. The JSON records the full command, digest parameters, reference denominator, sequencing budget, scoring weights, warnings, output paths, and ranked results for reproducibility.

The default `--objective balanced` prioritizes pairs that both match the requested weighted genome percentage within `--coverage-tolerance-pct` and meet `--desired-depth`/`--depth` for the planned sample count. Other objectives are available for sensitivity checks: `closest-coverage`, `depth-first`, `feasible-lowest-coverage`, and `max-depth`.

This inverse-design command is still a sequence-level model. It does not model methylation sensitivity, partial digestion, star activity, enzyme efficiency, buffer compatibility, empirical digestion rates, or per-locus depth dispersion.
