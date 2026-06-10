# radigest

[![CI](https://img.shields.io/github/actions/workflow/status/ericksamera/radigest/ci.yml?branch=main&label=ci)](https://github.com/ericksamera/radigest/actions/workflows/ci.yml)
[![DOI](http://zenodo.org/badge/979818941.svg)](https://doi.org/10.5281/zenodo.20176743)
[![Go](https://img.shields.io/badge/go-%3E=%201.22-blue)](https://golang.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)

`radigest` is a fast in-silico restriction digest and enzyme-pair screening toolkit for genomics.

It does two main things:

1. **Digest known enzymes quickly.**
   Use `radigest` when you already know the enzyme or enzyme pair.

2. **Screen enzyme pairs for a design target.**
   Use `radigest-design` when you know the genome fraction, sample count, read budget, and depth target, but still need to choose enzymes.

The model is deterministic and sequence-level. It finds recognition sites in a reference FASTA, applies digest rules, applies optional size-selection weights, and writes reproducible outputs.

---

## Install

```bash
make build
make install
```

This installs:

```text
radigest
radigest-design
radigest-fit-size-model
```

For development/helper commands:

```bash
make build-dev
make install-dev
```

---

## Which command should I use?

| Situation | Command |
|---|---|
| I already know my enzyme or enzyme pair | `radigest` |
| I want BED/GFF/TSV/FASTA fragment outputs | `radigest` |
| I want to screen enzyme pairs against a design target | `radigest-design` |
| I want to fit a size-selection model from observed inserts | `radigest-fit-size-model` |

---

# 1. Fast digestion with `radigest`

Use `radigest` when the enzyme choice is already known.

## Minimal digest

```bash
radigest -fasta ref.fa -enzymes EcoRI,MseI
```

With no output flags, `radigest` writes a JSON run summary to stdout.

Save the summary:

```bash
radigest -fasta ref.fa -enzymes EcoRI,MseI \
  -json run.json
```

## Size-select fragments

```bash
radigest \
  -fasta ref.fa \
  -enzymes EcoRI,MseI \
  -min 300 \
  -max 600 \
  -json run.json
```

## Write fragment files

```bash
radigest \
  -fasta ref.fa \
  -enzymes EcoRI,MseI \
  -min 300 \
  -max 600 \
  -bed fragments.bed \
  -gff fragments.gff3 \
  -fragments-tsv fragments.tsv \
  -json run.json
```

Common outputs:

| Output | Flag |
|---|---|
| JSON run summary | `-json run.json` |
| BED6 fragments | `-bed fragments.bed` |
| GFF3 fragments | `-gff fragments.gff3` |
| Per-fragment TSV | `-fragments-tsv fragments.tsv` |
| Fragment sequences | `-fragments-fasta fragments.fa` |

Coordinates:

| Format | Coordinates |
|---|---|
| GFF3 | 1-based closed |
| BED | 0-based half-open |
| TSV | 0-based half-open |
| FASTA metadata | 0-based half-open |

## Double-digest behavior

By default, double-digest mode keeps adjacent **AB/BA** fragments:

```bash
radigest -fasta ref.fa -enzymes EcoRI,MseI
```

To also keep **AA/BB** adjacent fragments:

```bash
radigest -fasta ref.fa -enzymes EcoRI,MseI -allow-same
```

Terminal contig-end fragments are omitted by default. Include them with:

```bash
radigest -fasta ref.fa -enzymes EcoRI,MseI -include-ends
```

## Size-selection models

The hard size window controls which fragments are retained:

```bash
-min 300 -max 600
```

For weighted recovery modeling, score a broader range:

```bash
radigest \
  -fasta ref.fa \
  -enzymes PstI,MspI \
  -min 300 \
  -max 600 \
  -score-min 1 \
  -score-max 2000 \
  -size-model normal \
  -size-mean 275 \
  -size-sd 85 \
  -fragments-tsv fragments.tsv \
  -json run.json
```

Supported models:

```text
hard
normal
triangular
soft-window
```

Use `hard` for a strict size window. Use the other models when size recovery is expected to be gradual rather than perfectly sharp.

---

# 2. Enzyme-pair screening with `radigest-design`

Use `radigest-design` when the experimental target is known and the enzyme pair is the unknown.

Typical question:

> Which enzyme pair best matches my target genome fraction and sequencing budget?

## Minimal design run

```bash
radigest-design \
  --ref ref.fa \
  --enzymes EcoRI,MseI,PstI,ApeKI,NlaIII,MspI \
  --pct 2.5 \
  --depth 10 \
  --samples 96 \
  --read-length 150 \
  --flowcell-read-pairs 300M \
  --usable-read-fraction 0.85
```

Aliases:

| Alias | Full flag |
|---|---|
| `--ref` | `--fasta` |
| `--pct` | `--target-genome-pct` |
| `--depth` | `--desired-depth` |

## Use an enzyme list

```bash
cat > candidate_enzymes.txt <<'EOF'
EcoRI
MseI
PstI
MspI
ApeKI
NlaIII
MluCI
BfaI
EOF
```

```bash
radigest-design \
  --ref ref.fa \
  --enzymes candidate_enzymes.txt \
  --pct 2.5 \
  --depth 10 \
  --samples 96 \
  --read-length 150 \
  --flowcell-read-pairs 300M \
  --usable-read-fraction 0.85 \
  --out-dir radigest_design
```

For broad exploration:

```bash
radigest-design ... --enzymes all
```

Review final enzyme choices manually before wet-lab use.

## Add size-selection assumptions

```bash
radigest-design \
  --ref ref.fa \
  --enzymes candidate_enzymes.txt \
  --pct 2.5 \
  --coverage-tolerance-pct 0.25 \
  --depth 10 \
  --samples 96 \
  --read-layout pe \
  --read-length 150 \
  --flowcell-read-pairs 300M \
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

## What `radigest-design` reports

`radigest-design` writes:

| File | Use |
|---|---|
| `design.summary.tsv` | Compact ranked table for human review |
| `design.tsv` | Full machine-readable ranked table |
| `design.json` | Full provenance and reproducibility record |
| `design.report.txt` | Simple key-value report |

It also prints a recommendation-first terminal summary:

```text
Recommendation:
Recommended pair: PstI,MspI
Status: feasible
Why: predicted 2.43% genome vs target 2.50%; predicted 12.1x mean locus depth vs target 10x
Budget: 96 samples, 300M read pairs, 0.85 usable fraction
Main caution: mean insert 247 bp is below 2x150 bp; paired-end overlap likely
Files: radigest_design/design.summary.tsv, radigest_design/design.report.txt
```

Start with:

```bash
column -ts $'\t' radigest_design/design.summary.tsv | less -S
```

## Key design terms

| Term | Meaning |
|---|---|
| `--pct` | Target weighted recovered genome percentage |
| `--depth` | Mean read-pair depth per recovered locus |
| `weighted_fragments` | Modeled recovered loci competing for reads |
| `predicted_depth` | Read pairs per sample divided by weighted fragments |
| `usable_read_fraction` | Fraction of read pairs expected to remain useful after demultiplexing/QC/deduplication |

`--depth` is not basewise WGS depth. It is a mean read-pair depth per recovered locus.

## Ranking objectives

Default:

```bash
--objective balanced
```

Other options:

```text
closest-coverage
depth-first
feasible-lowest-coverage
max-depth
```

Use `balanced` first. Rerun with another objective for sensitivity checks.

---

# Practical analysis workflow

## A. Digest a known pair

```bash
radigest \
  -fasta ref.fa \
  -enzymes PstI,MspI \
  -min 300 \
  -max 600 \
  -bed fragments.bed \
  -fragments-tsv fragments.tsv \
  -json run.json
```

## B. Choose a pair, then digest it

```bash
radigest-design \
  --ref ref.fa \
  --enzymes candidate_enzymes.txt \
  --pct 2.5 \
  --depth 10 \
  --samples 96 \
  --read-length 150 \
  --flowcell-read-pairs 300M \
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

Inspect the recommendation:

```bash
column -ts $'\t' radigest_design/design.summary.tsv | less -S
cat radigest_design/design.report.txt
```

Then digest the selected pair:

```bash
radigest \
  -fasta ref.fa \
  -enzymes PstI,MspI \
  -min 300 \
  -max 600 \
  -bed final_fragments.bed \
  -fragments-tsv final_fragments.tsv \
  -json final_digest.json
```

---

# Model scope

`radigest` models:

- recognition sites
- cut coordinates
- single- and double-digest fragment rules
- hard size windows
- optional size-selection weights
- weighted recovered genome percentage
- mean read-pair depth per recovered locus

It does **not** model:

- methylation sensitivity
- partial digestion
- star activity
- enzyme efficiency
- buffer compatibility
- empirical digestion rates
- per-locus depth dispersion

Enzymes with the same recognition motif and cut coordinate are treated identically by the sequence-level model.

Use `radigest` for fast, reproducible in-silico screening. Validate final enzyme choices against wet-lab constraints.

---

# Help

```bash
radigest --help
radigest-design --help
radigest-fit-size-model --help
radigest -list-enzymes
```

---

# Citation

If you use `radigest`, please cite the DOI listed at the top of this README.
