# radigest

[![CI](https://img.shields.io/github/actions/workflow/status/ericksamera/radigest/ci.yml?branch=main&label=ci)](https://github.com/ericksamera/radigest/actions/workflows/ci.yml)
[![DOI](http://zenodo.org/badge/979818941.svg)](https://doi.org/10.5281/zenodo.20176743)
[![Go](https://img.shields.io/badge/go-%3E=%201.22-blue)](https://golang.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)

Fast in-silico restriction digest for genomics. Give it a reference FASTA (plain or `.gz`) or synthesize one on the fly; pass one or two enzymes; it scans, size-selects, and exports fragments as GFF3 for GBS/ddRAD, probe design, or visualization. Output order is deterministic even with multithreading.

---

## Features

- **Single or double digest.** Double-digest keeps **adjacent AB/BA** by default; enable **AA/BB** too with `-allow-same`. Single-digest uses consecutive A cuts. Terminal chromosome/contig-end fragments are omitted by default; keep them with `-include-ends`.
- **IUPAC & cut offsets.** Sites accept degenerate codes; the cut index comes from `^` in the site (or mid-site if missing). `-strict-cuts` makes missing carets an error.
- **Robust FASTA I/O.** Read from a path or `-` (STDIN), auto-detect `.gz`, normalize case, and **trim CRLF**. `N` in the **reference** does **not** match any site.
- **Synthetic genomes.** Generate a single-chromosome genome named `chr1` with `-sim-len`, `-sim-gc`, `-sim-seed` and digest it directly—no FASTA on disk needed.
- **Clean outputs.** GFF3 with `ID=<chr>_<n>;Length=<bp>`; per-fragment TSV with insert length, hard-keep status, and size-selection weight; optional JSON summary of counts/bases per chromosome. Coordinates are **1-based closed** in GFF (internally 0-based half-open).
- **Size-selection scoring.** Keep the hard `-min/-max` window for GFF while assigning per-fragment recovery weights with `hard`, `normal`, `triangular`, or `soft-window` models over an optional broader score range.
- **Streaming fragment export.** The CLI streams digest fragments to the collector instead of materializing every kept fragment for a chromosome before writing.

---

## Quick start

```bash
# Single digest (EcoRI) → GFF file
radigest -fasta ref.fa -enzymes EcoRI -gff fragments.gff3

# Double digest with size selection + JSON summary
radigest -fasta ref.fa -enzymes EcoRI,MseI -min 100 -max 800 -gff fragments.gff3 -json run.json

# ddRAD-style soft-window scoring with broad per-fragment TSV for downstream modeling
radigest -fasta ref.fa -enzymes PstI,MspI -min 250 -max 500 -score-min 1 -score-max 1000 \
  -size-model soft-window -size-edge-sd 25 -fragments-tsv fragments.tsv -json run.json

# Double digest but ALSO keep AA/BB neighbors
radigest -fasta ref.fa -enzymes EcoRI,MseI -allow-same -gff fragments.gff3

# Include terminal fragments from chromosome/contig ends to the nearest cut
radigest -fasta ref.fa -enzymes EcoRI,MseI -include-ends -gff fragments.gff3

# Simulate a 10 Mb genome at 42% GC and digest
radigest -sim-len 10000000 -sim-gc 0.42 -sim-seed 123 -enzymes EcoRI,MseI -gff out.gff3
```

---

## CLI (most used)

- `-fasta <path|->` — reference FASTA; `-` = STDIN; `.gz` auto-detected.
- `-enzymes E1[,E2]` — one (A) or two (A,B) only. In double-digest, AB/BA by default.
- `-min/-max` — hard-selected insert-length window used for GFF output and `hard_kept` in TSV (**default min=1**).
- `-score-min/-score-max` — broader insert-length range emitted to `-fragments-tsv` and used for weighted size-selection stats; defaults to `-min/-max`.
- `-size-model hard|normal|triangular|soft-window` — per-fragment size-selection weight model (**default `hard`**).
- `-size-mean`, `-size-sd`, `-size-edge-sd` — parameters for `normal`, `triangular`, and `soft-window` scoring.
- `-gff <path|->` — GFF3 out for hard-kept fragments; `-` = STDOUT.
- `-fragments-tsv <path|->` — per-fragment TSV for score-range fragments (**default `fragments.tsv`; empty string disables**).
- `-json <path>` — write a run summary (counts, bases, per-chrom stats, and size-selection weighted stats).
- `-threads <n>` — positive worker count; `-v`, `-version`, `-list-enzymes`.
- **Simulation:** `-sim-len <bp>`, `-sim-gc <0..1>` (invalid values error), `-sim-seed <int>` (emits a single `chr1`).
- **Modes:** `-allow-same` (keep AA/BB in double-digest), `-include-ends` (also emit terminal chromosome/contig-end fragments), `-strict-cuts` (error if a site lacks `^` and would otherwise fall back to mid-site).

---

## Scope and limitations

radigest is a deterministic sequence-level model. It identifies recognition sites and cut coordinates from the reference sequence only. It does **not** model methylation sensitivity, partial digestion, star activity, enzyme efficiency, buffer compatibility, or empirical digestion rates. Enzymes with the same recognition motif and cut coordinate are treated identically by the digest logic even when their wet-lab behavior can differ under methylation or assay conditions.

---

## Outputs

### GFF3

```
##gff-version 3
chr1	radigest	fragment	<start>	<end>	.	+	.	ID=chr1_1;Length=123
```

`start/end` are **1-based closed**; `Length` is `end - start + 1`. Ordering is deterministic per chromosome.

### Fragment TSV

By default, radigest also writes a per-fragment TSV for fragments in the score range:

```
chrom	start0	end0	length	hard_kept	size_weight
chr1	10422	10731	309	true	0.982143
chr1	18831	18922	91	false	0.014221
```

`hard_kept` is true when the insert length is inside `-min/-max`. `size_weight` is the selected size model evaluated on insert length only. Use `-fragments-tsv ""` to disable this output.

### JSON (optional)

```json
{
  "enzymes": ["EcoRI", "MseI"],
  "min_length": 100,
  "max_length": 800,
  "fragments_tsv": "fragments.tsv",
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

The screen writes one JSON summary per pair under `pair_screen/json/` and logs under `pair_screen/logs/`. It disables GFF and fragment TSV output for speed during initial screening.

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
```

The ranked TSV includes `weighted_bases`, `weighted_fragments`, `raw_bases_in_window`, `raw_fragments_in_window`, `mean_weighted_length`, and genome-percentage columns when a denominator is supplied.
