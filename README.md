# radigest

Fast in‑silico restriction digest for genomics. Give it a reference FASTA (`.fna(.gz)`) and one or two enzyme names; it scans the genome, applies size selection, and exports fragments as GFF3 for downstream workflows (GBS/ddRAD, probe design, visualization). Fragments are produced deterministically even with multithreading.&#x20;

---

## Features

* **Single or double digest.** In double‑digest mode only **adjacent AB/BA** fragments are kept; single‑digest uses consecutive A cuts.&#x20;
* **IUPAC recognition & cut offsets.** Enzyme sites accept degenerate bases; cut position comes from `^` in the site or defaults to mid‑site.&#x20;
* **Stream‑friendly + fast.** Reads FASTA from a file or `-` (STDIN), auto‑detects `.gz`, uses a worker pool, and writes fragments in a stable order.&#x20;
* **Clean outputs.** GFF3 lines with `ID=<chr>_<n>;Length=<bp>` attributes; optional machine‑readable JSON summary. Coordinates in GFF are **1‑based closed**.&#x20;
* **Practical CLI.** `-min/-max` for size selection; `-list-enzymes` to discover supported names; `-threads` and `-v` for control.&#x20;

---

## Quick start

```bash
# Single digest (EcoRI) → GFF file
radigest -fasta ref.fa -enzymes EcoRI -gff fragments.gff3

# Double digest (EcoRI + MseI), size‑select 100–800 bp, and write a JSON summary
radigest -fasta ref.fa -enzymes EcoRI,MseI -min 100 -max 800 -gff fragments.gff3 -json run.json

# Stream a compressed FASTA in, write GFF to stdout
zcat ref.fa.gz | radigest -fasta - -enzymes EcoRI,MseI -gff -
```

Enzyme names come from the built‑in database (e.g., `EcoRI`, `MseI`, `MspI`, `HindIII`). Use `-list-enzymes` to see them. The **first two names** define the AB pair; they must differ. Additional names, if provided, are ignored today.&#x20;

---

## Command‑line options (most used)

* `-fasta <path|->` — reference FASTA; `-` means STDIN; `.gz` auto‑detected.&#x20;
* `-enzymes E1[,E2]` — one (A) or two (A,B) enzymes; AB/BA adjacency enforced in double‑digest.&#x20;
* `-min/-max` — keep fragments in `[min,max]` bp.&#x20;
* `-gff <path|->` — GFF3 output file or `-` for STDOUT.&#x20;
* `-json <path>` — optional run summary (see below).&#x20;
* `-threads <n>`, `-v`, `-version`, `-list-enzymes`.&#x20;

---

## Outputs

### GFF3

One header plus one line per kept fragment:

```
##gff-version 3
chr1  radigest  fragment  <start>  <end>  .  +  .  ID=chr1_1;Length=123
```

`start/end` are **1‑based closed**; `Length` is the bp span. Ordering is deterministic per chromosome.&#x20;

### JSON summary

Written when `-json` is set:

```json
{
  "enzymes": ["EcoRI","MseI"],
  "min_length": 100,
  "max_length": 800,
  "total_fragments": 123456,
  "total_bases": 7891011,
  "per_chromosome": {
    "chr1": {"fragments": 23456, "bases": 3456789}
  }
}
```

Field names match the program’s embedded stats (`total_fragments`, `total_bases`, `per_chromosome[chr].{fragments,bases}`).&#x20;

---

## Notes for bioinformatics use

* Intended for **GBS/ddRAD** style in‑silico selection and probe/amplicon panel design.
* Internally uses 0‑based half‑open positions; GFF is 1‑based closed.&#x20;
* IUPAC matching enables methylation‑insensitive patterns and degenerate sites; cut offsets honor `^`.&#x20;
