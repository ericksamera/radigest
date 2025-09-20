# radigest

Fast in-silico restriction digest for genomics. Give it a reference FASTA (plain or `.gz`) or synthesize one on the fly; pass one or two enzymes; it scans, size-selects, and exports fragments as GFF3 for GBS/ddRAD, probe design, or visualization. Output order is deterministic even with multithreading.

---

## Features

* **Single or double digest.** Double-digest keeps **adjacent AB/BA** by default; enable **AA/BB** too with `-allow-same`. Single-digest uses consecutive A cuts.
* **IUPAC & cut offsets.** Sites accept degenerate codes; the cut index comes from `^` in the site (or mid-site if missing). `-strict-cuts` makes missing carets an error.
* **Robust FASTA I/O.** Read from a path or `-` (STDIN), auto-detect `.gz`, normalize case, and **trim CRLF**. `N` in the **reference** does **not** match any site.
* **Synthetic genomes.** Generate a single-chromosome genome named `chr1` with `-sim-len`, `-sim-gc`, `-sim-seed` and digest it directly—no FASTA on disk needed.
* **Clean outputs.** GFF3 with `ID=<chr>_<n>;Length=<bp>`; optional JSON summary of counts/bases per chromosome. Coordinates are **1-based closed** in GFF (internally 0-based half-open).

---

## Quick start

```bash
# Single digest (EcoRI) → GFF file
radigest -fasta ref.fa -enzymes EcoRI -gff fragments.gff3

# Double digest with size selection + JSON summary
radigest -fasta ref.fa -enzymes EcoRI,MseI -min 100 -max 800 -gff fragments.gff3 -json run.json

# Double digest but ALSO keep AA/BB neighbors
radigest -fasta ref.fa -enzymes EcoRI,MseI -allow-same -gff fragments.gff3

# Simulate a 10 Mb genome at 42% GC and digest
radigest -sim-len 10000000 -sim-gc 0.42 -sim-seed 123 -enzymes EcoRI,MseI -gff out.gff3
```

---

## CLI (most used)

* `-fasta <path|->` — reference FASTA; `-` = STDIN; `.gz` auto-detected.
* `-enzymes E1[,E2]` — one (A) or two (A,B). In double-digest, AB/BA by default.
* `-min/-max` — keep fragments in `[min,max]` bp (**default min=1**).
* `-gff <path|->` — GFF3 out; `-` = STDOUT.
* `-json <path>` — write a run summary (counts, bases, per-chrom stats).
* `-threads <n>`, `-v`, `-version`, `-list-enzymes`.
* **Simulation:** `-sim-len <bp>`, `-sim-gc <0..1>`, `-sim-seed <int>` (emits a single `chr1`).
* **Modes:** `-allow-same` (keep AA/BB in double-digest), `-strict-cuts` (error if a site lacks `^` and would otherwise fall back to mid-site).

---

## Outputs

### GFF3

```
##gff-version 3
chr1	radigest	fragment	<start>	<end>	.	+	.	ID=chr1_1;Length=123
```

`start/end` are **1-based closed**; `Length` is `end - start + 1`. Ordering is deterministic per chromosome.

### JSON (optional)

```json
{
  "enzymes": ["EcoRI","MseI"],
  "min_length": 100,
  "max_length": 800,
  "total_fragments": 123456,
  "total_bases": 7891011,
  "per_chromosome": {"chr1": {"fragments": 23456, "bases": 3456789}}
}
```