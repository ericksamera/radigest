package main

import (
	"fmt"
	"io"
	"runtime"
	"strconv"

	"github.com/ericksamera/radigest/internal/clihelp"
	"github.com/ericksamera/radigest/internal/design"
)

func writeDesignUsage(w io.Writer, version string, weights design.ScoreWeights) {
	_, _ = fmt.Fprintln(w, "radigest-design")
	_, _ = fmt.Fprintln(w, "Author:  Gennerick J. Samera (erick.samera@kpu.ca)")
	_, _ = fmt.Fprintf(w, "Version: %s\n", version)
	_, _ = fmt.Fprintln(w, "License: MIT")
	_, _ = fmt.Fprintln(w, "Description: rank enzyme pairs from genome-fraction, depth, and sequencing-budget targets")
	_, _ = fmt.Fprintln(w)
	_, _ = fmt.Fprintln(w, "Usage:")
	_, _ = fmt.Fprintln(w, "  radigest-design --ref ref.fa --enzymes EcoRI,MseI,PstI,ApeKI --pct 2.5 --depth 10 --samples 96 --read-length 150 --flowcell-read-pairs 300M [options]")
	_, _ = fmt.Fprintln(w, "  radigest-design --fasta ref.fa --enzymes candidate_enzymes.txt --target-genome-pct 2.5 --desired-depth 10 --samples 96 --read-length 150 --lane-read-pairs 300M --lanes 1 [options]")
	_, _ = fmt.Fprintln(w)

	clihelp.WriteFlagGroups(w, []clihelp.Group{
		{
			Title: "Required design inputs",
			Items: []clihelp.Flag{
				{Names: []string{"--ref", "--fasta"}, Arg: "PATH", Text: "Reference FASTA. Plain or .gz."},
				{Names: []string{"--enzymes"}, Arg: "LIST|FILE|all", Text: "Candidate enzymes as comma-separated names, a one-per-line file, or 'all'."},
				{Names: []string{"--pct", "--target-genome-pct"}, Arg: "FLOAT", Text: "Target weighted genome percentage, for example 2.5."},
				{Names: []string{"--depth", "--target-depth", "--desired-depth"}, Arg: "FLOAT", Text: "Target mean read-pair depth per recovered locus."},
				{Names: []string{"--samples"}, Arg: "INT", Text: "Planned number of samples."},
				{Names: []string{"--read-length"}, Arg: "INT", Text: "Sequencing read length in bp."},
				{Names: []string{"--flowcell-read-pairs"}, Arg: "COUNT", Text: "Total read pairs for a one-flowcell/run budget, for example 50M or 300M. Mutually exclusive with --lane-read-pairs."},
				{Names: []string{"--lane-read-pairs"}, Arg: "COUNT", Text: "Read pairs per lane, for example 300M. Mutually exclusive with --flowcell-read-pairs."},
			},
		},
		{
			Title: "Size selection and recovery model",
			Items: []clihelp.Flag{
				{Names: []string{"--min"}, Arg: "INT", Default: "300", Text: "Hard lower insert-size bound in bp."},
				{Names: []string{"--max"}, Arg: "INT", Default: "600", Text: "Hard upper insert-size bound in bp."},
				{Names: []string{"--score-min"}, Arg: "INT", Default: "1", Text: "Lower insert-size bound included in recovery-weight scoring."},
				{Names: []string{"--score-max"}, Arg: "INT", Default: "2000", Text: "Upper insert-size bound included in recovery-weight scoring."},
				{Names: []string{"--size-model"}, Arg: "MODEL", Default: "normal", Text: "Size-selection/recovery weighting model: hard, normal, triangular, or soft-window."},
				{Names: []string{"--size-mean"}, Arg: "FLOAT", Default: "275", Text: "Peak/target insert length for normal and triangular models."},
				{Names: []string{"--size-sd"}, Arg: "FLOAT", Default: "85", Text: "Standard deviation for the normal model."},
				{Names: []string{"--size-edge-sd"}, Arg: "FLOAT", Default: "25", Text: "Edge softness for the soft-window model."},
			},
		},
		{
			Title: "Sequencing budget",
			Items: []clihelp.Flag{
				{Names: []string{"--lanes"}, Arg: "INT", Default: "1", Text: "Number of lanes. Only valid with --lane-read-pairs."},
				{Names: []string{"--usable-read-fraction"}, Arg: "FLOAT", Default: "1", Text: "Fraction of read pairs usable after demultiplexing/QC/deduplication."},
				{Names: []string{"--read-layout"}, Arg: "pe|se", Default: "pe", Text: "Read layout used for insert diagnostics."},
			},
		},
		{
			Title: "Reference denominator",
			Items: []clihelp.Flag{
				{Names: []string{"--denominator"}, Arg: "non-n|all", Default: "non-n", Text: "Denominator used for genome-percentage calculations."},
				{Names: []string{"--genome-bases"}, Arg: "COUNT", Text: "Explicit denominator. Skips FASTA base counting."},
			},
		},
		{
			Title: "Digest behavior",
			Items: []clihelp.Flag{
				{Names: []string{"--allow-same"}, Text: "In double-digest scoring, also keep AA/BB adjacent fragments."},
				{Names: []string{"--include-ends"}, Text: "Include terminal fragments from contig ends to nearest cut."},
				{Names: []string{"--strict-cuts"}, Text: "Error if an enzyme lacks an explicit cut coordinate."},
			},
		},
		{
			Title: "Ranking and scoring",
			Items: []clihelp.Flag{
				{Names: []string{"--coverage-tolerance-pct"}, Arg: "FLOAT", Default: "0.25", Text: "Absolute tolerance around --target-genome-pct."},
				{Names: []string{"--objective"}, Arg: "OBJECTIVE", Default: "balanced", Text: "Ranking objective: balanced, closest-coverage, depth-first, feasible-lowest-coverage, or max-depth."},
				{Names: []string{"--weight-coverage"}, Arg: "FLOAT", Default: formatHelpFloat(weights.Coverage), Text: "Fit-loss weight for coverage error."},
				{Names: []string{"--weight-depth"}, Arg: "FLOAT", Default: formatHelpFloat(weights.Depth), Text: "Fit-loss weight for depth shortfall."},
				{Names: []string{"--weight-overcoverage"}, Arg: "FLOAT", Default: formatHelpFloat(weights.Overcoverage), Text: "Additional fit-loss weight for overcoverage."},
				{Names: []string{"--weight-insert"}, Arg: "FLOAT", Default: formatHelpFloat(weights.Insert), Text: "Fit-loss weight for insert-size risk."},
			},
		},
		{
			Title: "Outputs",
			Items: []clihelp.Flag{
				{Names: []string{"--out-dir"}, Arg: "DIR", Default: "radigest_design", Text: "Output directory."},
				{Names: []string{"--summary-tsv"}, Arg: "PATH", Default: "<out-dir>/design.summary.tsv", Text: "Compact human-review table."},
				{Names: []string{"--tsv"}, Arg: "PATH", Default: "<out-dir>/design.tsv", Text: "Full machine-readable table."},
				{Names: []string{"--json"}, Arg: "PATH", Default: "<out-dir>/design.json", Text: "Reproducibility/report JSON."},
				{Names: []string{"--top"}, Arg: "INT", Default: "0", Text: "Limit reported ranked rows/results. 0 means all."},
				{Names: []string{"--force"}, Text: "Overwrite existing outputs."},
			},
		},
		{
			Title: "Performance",
			Items: []clihelp.Flag{
				{Names: []string{"--jobs"}, Arg: "INT", Default: "--threads", Text: "Parallel pair-scoring workers."},
				{Names: []string{"--threads"}, Arg: "INT", Default: fmt.Sprintf("%d", runtime.NumCPU()), Text: "Worker-count fallback when --jobs is unset."},
				{Names: []string{"--build-workers"}, Arg: "INT", Default: "--jobs, then --threads", Text: "Parallel cut-index build workers."},
				{Names: []string{"--max-pairs"}, Arg: "INT", Default: "0", Text: "Score at most this many candidate pairs. 0 means all."},
			},
		},
		{
			Title: "Other",
			Items: []clihelp.Flag{
				{Names: []string{"--version"}, Text: "Print version and exit."},
				{Names: []string{"--help", "-h"}, Text: "Show this help."},
			},
		},
	})

	clihelp.WriteSizeModelReference(w, "--")

	_, _ = fmt.Fprintln(w, "Outputs written:")
	_, _ = fmt.Fprintln(w, "  radigest_design/design.summary.tsv  compact table for human review")
	_, _ = fmt.Fprintln(w, "  radigest_design/design.tsv          full ranked design table")
	_, _ = fmt.Fprintln(w, "  radigest_design/design.json         full provenance, parameters, warnings, and ranked results")
	_, _ = fmt.Fprintln(w)
	_, _ = fmt.Fprintln(w, "Notes:")
	_, _ = fmt.Fprintln(w, "  Genome percentage means weighted recovered genome percentage under the specified size-selection/recovery model.")
	_, _ = fmt.Fprintln(w, "  Depth means mean read-pair depth per recovered locus, not basewise WGS depth.")
	_, _ = fmt.Fprintln(w, "  The model is sequence-level only; it does not model methylation sensitivity, partial digestion, enzyme efficiency, buffer compatibility, or per-locus depth dispersion.")
}

func formatHelpFloat(v float64) string {
	return strconv.FormatFloat(v, 'f', -1, 64)
}
