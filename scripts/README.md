# radigest helper commands

These pure-Python helpers are intended to be installed alongside the Go
`radigest` binary, for example by Bioconda:

```bash
radigest --help
radigest-screen-pairs --help
radigest-rank-pairs --help
radigest-fit-size-model --help
```

## Commands

### `radigest-screen-pairs`

Runs `radigest` for all unique enzyme pairs from a candidate list and writes one
JSON summary per pair.

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
