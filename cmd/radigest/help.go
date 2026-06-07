package main

import (
	"fmt"
	"io"
	"runtime"

	"github.com/ericksamera/radigest/internal/clihelp"
)

func writeRadigestUsage(w io.Writer, version string) {
	_, _ = fmt.Fprintln(w, "radigest")
	_, _ = fmt.Fprintln(w, "Author:  Gennerick J. Samera (erick.samera@kpu.ca)")
	_, _ = fmt.Fprintf(w, "Version: %s\n", version)
	_, _ = fmt.Fprintln(w, "License: MIT")
	_, _ = fmt.Fprintln(w, "Description: deterministic in-silico restriction digest")
	_, _ = fmt.Fprintln(w)
	_, _ = fmt.Fprintln(w, "Usage:")
	_, _ = fmt.Fprintln(w, "  radigest -fasta <ref.fa|-> -enzymes <E1[,E2]> [options]")
	_, _ = fmt.Fprintln(w, "  radigest -sim-len <bp> -sim-gc <0..1> -enzymes <E1[,E2]> [options]")
	_, _ = fmt.Fprintln(w)

	clihelp.WriteFlagGroups(w, []clihelp.Group{
		{
			Title: "Required inputs",
			Intro: []string{"Provide -enzymes and exactly one of -fasta or -sim-len."},
			Items: []clihelp.Flag{
				{Names: []string{"-enzymes"}, Arg: "E1[,E2]", Text: "One or two enzyme names. Single digest uses consecutive A cuts. Double digest keeps adjacent AB/BA fragments by default."},
				{Names: []string{"-fasta"}, Arg: "PATH|-", Text: "Reference FASTA. Plain, .gz, or '-' for stdin."},
				{Names: []string{"-sim-len"}, Arg: "BP", Text: "Simulate a single chromosome named chr1 instead of reading FASTA."},
			},
		},
		{
			Title: "Digest behavior",
			Items: []clihelp.Flag{
				{Names: []string{"-min"}, Arg: "INT", Default: "1", Text: "Hard lower insert-size bound for retained fragments."},
				{Names: []string{"-max"}, Arg: "INT", Default: "1073741824", Text: "Hard upper insert-size bound for retained fragments."},
				{Names: []string{"-allow-same"}, Text: "In double digest, also keep AA/BB adjacent fragments."},
				{Names: []string{"-include-ends"}, Text: "Include terminal fragments from contig ends to the nearest cut."},
				{Names: []string{"-strict-cuts"}, Text: "Error if an enzyme lacks an explicit cut coordinate."},
			},
		},
		{
			Title: "Size filtering and scoring",
			Items: []clihelp.Flag{
				{Names: []string{"-score-min"}, Arg: "INT", Default: "-min", Text: "Lower insert-size bound included in size-selection scoring and fragment TSV output."},
				{Names: []string{"-score-max"}, Arg: "INT", Default: "-max", Text: "Upper insert-size bound included in size-selection scoring and fragment TSV output."},
				{Names: []string{"-size-model"}, Arg: "MODEL", Default: "hard", Text: "Size-selection/recovery weighting model: hard, normal, triangular, or soft-window."},
				{Names: []string{"-size-mean"}, Arg: "FLOAT", Default: "midpoint of -min/-max", Text: "Peak/target insert length for normal and triangular models."},
				{Names: []string{"-size-sd"}, Arg: "FLOAT", Default: "35", Text: "Standard deviation for the normal model."},
				{Names: []string{"-size-edge-sd"}, Arg: "FLOAT", Default: "25", Text: "Edge softness for the soft-window model."},
			},
		},
		{
			Title: "Outputs",
			Intro: []string{"If no output flags are set, run summary JSON is written to stdout."},
			Items: []clihelp.Flag{
				{Names: []string{"-json"}, Arg: "PATH|-", Text: "Run summary JSON."},
				{Names: []string{"-gff"}, Arg: "PATH|-", Text: "GFF3 for hard-kept fragments."},
				{Names: []string{"-bed"}, Arg: "PATH|-", Text: "BED6 for hard-kept fragments."},
				{Names: []string{"-fragments-tsv"}, Arg: "PATH|-", Text: "Per-fragment TSV for score-range fragments."},
				{Names: []string{"-fragments-fasta"}, Arg: "PATH|-", Text: "FASTA sequences for hard-kept fragments."},
			},
		},
		{
			Title: "Performance",
			Items: []clihelp.Flag{
				{Names: []string{"-threads"}, Arg: "INT", Default: fmt.Sprintf("%d", runtime.NumCPU()), Text: "Number of worker goroutines."},
				{Names: []string{"-v"}, Text: "Print verbose progress to stderr."},
			},
		},
		{
			Title: "Simulation",
			Items: []clihelp.Flag{
				{Names: []string{"-sim-gc"}, Arg: "FLOAT", Default: "0.50", Text: "Target GC fraction in [0,1] for -sim-len."},
				{Names: []string{"-sim-seed"}, Arg: "INT", Default: "1", Text: "PRNG seed for -sim-len. Use 0 for a time-based seed; the resolved seed is recorded in JSON."},
			},
		},
		{
			Title: "Other",
			Items: []clihelp.Flag{
				{Names: []string{"-list-enzymes"}, Text: "List available enzyme names and exit."},
				{Names: []string{"-version"}, Text: "Print version and exit."},
				{Names: []string{"-help", "-h"}, Text: "Show this help."},
			},
		},
	})

	clihelp.WriteSizeModelReference(w, "-")

	_, _ = fmt.Fprintln(w, "Examples:")
	_, _ = fmt.Fprintln(w, "  radigest -fasta ref.fa -enzymes EcoRI,MseI")
	_, _ = fmt.Fprintln(w, "  zcat ref.fa.gz | radigest -fasta - -enzymes EcoRI,MseI -bed - | bgzip > frag.bed.gz")
	_, _ = fmt.Fprintln(w, "  radigest -fasta ref.fa -enzymes PstI,MspI -min 250 -max 500 -score-min 1 -score-max 1000 -size-model soft-window -size-edge-sd 25 -fragments-tsv fragments.tsv -json run.json")
	_, _ = fmt.Fprintln(w, "  radigest -sim-len 10000000 -sim-gc 0.42 -enzymes EcoRI,MseI -gff out.gff3")
	_, _ = fmt.Fprintln(w)
	_, _ = fmt.Fprintln(w, "Notes:")
	_, _ = fmt.Fprintln(w, "  Coordinates are 1-based closed in GFF and 0-based half-open in BED, TSV, and FASTA metadata.")
	_, _ = fmt.Fprintln(w, "  The model is sequence-level only; it does not model methylation sensitivity, partial digestion, enzyme efficiency, or buffer compatibility.")
}
