# radigest helper commands

These pure-Python helpers are intended to be installed alongside the Go
`radigest` binary, for example by Bioconda:

```bash
radigest --help
radigest-screen-pairs --help
radigest-rank-pairs --help
radigest-plan-depth --help
radigest-fit-size-model --help
```

## Commands

### `radigest-screen-pairs`

Runs `radigest` for all unique enzyme pairs from a candidate list and writes one
JSON summary per pair. GFF, TSV, and FASTA artifact outputs are not requested
during screening.

```bash
radigest-screen-pairs \
  --fasta ref.fa \
  --enzymes candidate_enzymes.txt \
  --min 300 \
  --max 600 \
  --score-min 1 \
  --score-max 2000 \
  --size-model normal \
  --size-mean 275 \
  --size-sd 85 \
  --out-dir pair_screen
```

When running from a source checkout instead of an installed package, pass:

```bash
--radigest 'go run ./cmd/radigest'
```

### `radigest-rank-pairs`

Ranks the JSON summaries produced by `radigest-screen-pairs`.

```bash
radigest-rank-pairs 'pair_screen/json/*.json' \
  --fasta ref.fa \
  --objective weighted-genome-pct \
  --out pair_screen/ranked_pairs.tsv

# Reuse a known non-N genome denominator without rereading the FASTA
GENOME_BASES=2643888753
radigest-rank-pairs 'pair_screen/json/*.json' \
  --genome-bases "$GENOME_BASES" \
  --objective closest-target \
  --target-genome-pct 1.5 \
  --out pair_screen/ranked_pairs.closest_1.5pct.tsv
```

### `radigest-plan-depth`

Appends deterministic sequencing-budget estimates to the TSV produced by
`radigest-rank-pairs`. This first-pass planner uses the size-weighted fragment
count as the number of recovered loci competing for reads. `--desired-depth` is
a mean read-pair depth per recovered locus, not a callable-depth distribution
model or basewise WGS depth. It reports expected mean depth for a planned sample
count, the number of samples supportable at a desired depth, and the
budget-supported weighted genome percentage. It does not model per-locus depth
dispersion.

```bash
radigest-plan-depth pair_screen/ranked_pairs.genome_pct.tsv \
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
radigest-plan-depth pair_screen/ranked_pairs.genome_pct.tsv \
  --read-layout pe \
  --read-length 150 \
  --flowcell-read-pairs 50M \
  --usable-read-fraction 0.85 \
  --samples 37 \
  --desired-depth 30 \
  --target-genome-pct 1.5 \
  --out pair_screen/depth_plan.flowcell.tsv
```

### `radigest-fit-size-model`

Fits simple empirical insert-length recovery models from a radigest fragments TSV
and one observed TLEN per line.

```bash
radigest-fit-size-model \
  --fragments fragments.tsv \
  --tlens all.tlen.tsv \
  --min 300 \
  --max 600 \
  --score-min 1 \
  --score-max 2000 \
  --out size_model_rankings.tsv
```

The fitted model is an empirical recovery curve. It includes effects from size
selection, short-fragment representation, PCR, sequencing, and mapping, so do
not interpret it as a pure wet-lab size-selection probability.
