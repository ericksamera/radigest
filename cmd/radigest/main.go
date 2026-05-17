package main

import (
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"io"
	"os"
	"runtime"
	"sort"
	"strings"
	"sync"

	"github.com/ericksamera/radigest/internal/collector"
	"github.com/ericksamera/radigest/internal/digest"
	"github.com/ericksamera/radigest/internal/enzyme"
	"github.com/ericksamera/radigest/internal/fasta"
	"github.com/ericksamera/radigest/internal/fragmentfasta"
	"github.com/ericksamera/radigest/internal/fragmenttsv"
	"github.com/ericksamera/radigest/internal/sim"
	"github.com/ericksamera/radigest/internal/sizeselect"
)

var (
	version = "v0.1.0"
)

const summarySchemaVersion = 1

type digestResult struct {
	idx    int
	chr    string
	seq    []byte
	frags  <-chan digest.Fragment
	errors <-chan error
}

type runSummary struct {
	SchemaVersion   int              `json:"schema_version"`
	RadigestVersion string           `json:"radigest_version"`
	Command         []string         `json:"command"`
	Input           inputSummary     `json:"input"`
	Parameters      parameterSummary `json:"parameters"`
	Outputs         outputSummary    `json:"outputs"`
	Warnings        []string         `json:"warnings"`

	// Backward-compatible top-level fields retained for existing downstream tools.
	Enzymes        []string         `json:"enzymes"`
	MinLength      int              `json:"min_length"`
	MaxLength      int              `json:"max_length"`
	GFF            string           `json:"gff,omitempty"`
	FragmentsTSV   string           `json:"fragments_tsv,omitempty"`
	FragmentsFASTA string           `json:"fragments_fasta,omitempty"`
	SizeSelection  sizeselect.Stats `json:"size_selection"`
	collector.Stats
}

type inputSummary struct {
	Source           string   `json:"source"`
	FASTA            string   `json:"fasta,omitempty"`
	SimLength        int      `json:"sim_length,omitempty"`
	SimGC            *float64 `json:"sim_gc,omitempty"`
	SimSeedRequested *int64   `json:"sim_seed_requested,omitempty"`
	SimSeedResolved  *int64   `json:"sim_seed_resolved,omitempty"`
}

type parameterSummary struct {
	MinLength   int     `json:"min_length"`
	MaxLength   int     `json:"max_length"`
	ScoreMin    int     `json:"score_min"`
	ScoreMax    int     `json:"score_max"`
	SizeModel   string  `json:"size_model"`
	SizeMean    float64 `json:"size_mean,omitempty"`
	SizeSD      float64 `json:"size_sd,omitempty"`
	SizeEdgeSD  float64 `json:"size_edge_sd,omitempty"`
	Threads     int     `json:"threads"`
	AllowSame   bool    `json:"allow_same"`
	StrictCuts  bool    `json:"strict_cuts"`
	IncludeEnds bool    `json:"include_ends"`
}

type outputSummary struct {
	JSON           string `json:"json,omitempty"`
	GFF            string `json:"gff,omitempty"`
	FragmentsTSV   string `json:"fragments_tsv,omitempty"`
	FragmentsFASTA string `json:"fragments_fasta,omitempty"`
}

type usageError struct {
	err error
}

func (e usageError) Error() string {
	return e.err.Error()
}

func (e usageError) Unwrap() error {
	return e.err
}

func main() {
	if err := run(os.Args[1:], os.Stdin, os.Stdout, os.Stderr); err != nil {
		fmt.Fprintln(os.Stderr, "error:", err)
		os.Exit(exitCode(err))
	}
}

func exitCode(err error) int {
	var usage usageError
	if errors.As(err, &usage) {
		return 2
	}
	return 1
}

func run(args []string, stdin io.Reader, stdout, stderr io.Writer) error {
	if stdin == nil {
		stdin = strings.NewReader("")
	}
	if stdout == nil {
		stdout = io.Discard
	}
	if stderr == nil {
		stderr = io.Discard
	}

	fs := flag.NewFlagSet("radigest", flag.ContinueOnError)
	fs.SetOutput(stderr)

	// ---- CLI flags ----------------------------------------------------------
	fastaPath := fs.String("fasta", "", "reference FASTA file")
	enzFlag := fs.String("enzymes", "", "comma-separated enzyme names (one or two; two form the AB pair)")
	minLen := fs.Int("min", 1, "minimum fragment length (bp) for hard-selected outputs")
	maxLen := fs.Int("max", 1<<30, "maximum fragment length (bp) for hard-selected outputs")
	gffPath := fs.String("gff", "", "optional GFF3 output for hard-selected fragments (path or '-' for stdout); empty string disables")
	fragmentsTSVPath := fs.String("fragments-tsv", "", "optional per-fragment TSV for score-range fragments (path or '-' for stdout); empty string disables")
	fragmentsFASTAPath := fs.String("fragments-fasta", "", "optional FASTA output for hard-selected fragments (path or '-' for stdout); empty string disables")
	jsonPath := fs.String("json", "", "optional run summary JSON output (path or '-' for stdout); if no output flags are set, JSON is written to stdout")
	threads := fs.Int("threads", runtime.NumCPU(), "number of worker goroutines")
	verbose := fs.Bool("v", false, "verbose progress to stderr")
	showVer := fs.Bool("version", false, "print version and exit")
	listEns := fs.Bool("list-enzymes", false, "list available enzyme names and exit")

	// size-selection scoring
	scoreMinFlag := fs.Int("score-min", -1, "minimum fragment length included in fragments TSV and size-selection stats; default -min")
	scoreMaxFlag := fs.Int("score-max", -1, "maximum fragment length included in fragments TSV and size-selection stats; default -max")
	sizeModel := fs.String("size-model", "hard", "size-selection model: hard, normal, triangular, or soft-window")
	sizeMean := fs.Float64("size-mean", 0, "target/peak insert length for normal/triangular models; default midpoint of -min/-max")
	sizeSD := fs.Float64("size-sd", 35, "standard deviation for -size-model normal")
	sizeEdgeSD := fs.Float64("size-edge-sd", 25, "edge softness for -size-model soft-window")

	// digest behavior & validation
	allowSame := fs.Bool("allow-same", false, "double digest: also keep AA/BB neighbors (default AB/BA only)")
	includeEnds := fs.Bool("include-ends", false, "also emit terminal fragments from chromosome/contig ends to the nearest cut")
	strictCuts := fs.Bool("strict-cuts", false, "error if an enzyme lacks a caret and CutIndex==0 (no mid-site fallback)")

	// synthetic genome flags
	simLen := fs.Int("sim-len", 0, "synthesize a single-chromosome genome of this length (bp) instead of reading -fasta")
	simGC := fs.Float64("sim-gc", 0.50, "target GC fraction in [0,1] for -sim-len")
	simSeed := fs.Int64("sim-seed", 1, "PRNG seed for -sim-len (0 ⇒ time-based)")

	fs.Usage = func() {
		b := &strings.Builder{}
		fmt.Fprintln(b, "Author:  Erick Samera (erick.samera@kpu.ca)")
		fmt.Fprintln(b, "License: MIT")
		fmt.Fprintln(b, "Version:", version)
		fmt.Fprintln(b)
		fmt.Fprintln(b, "radigest — in-silico single/double digest with JSON summaries and optional fragment exports")
		fmt.Fprintln(b)
		fmt.Fprintln(b, "Usage:")
		fmt.Fprintln(b, "  radigest -fasta <ref.fa|-> -enzymes <E1[,E2]> [options]")
		fmt.Fprintln(b, "  radigest -sim-len <bp> -sim-gc <0..1> -enzymes <E1[,E2]> [options]")
		fmt.Fprintln(b)
		fmt.Fprintln(b, "Required flags:")
		fmt.Fprintln(b, "  -enzymes, and exactly one of -fasta or -sim-len")
		fmt.Fprintln(b)
		fmt.Fprintln(b, "Options:")
		oldOutput := fs.Output()
		fs.SetOutput(b)
		fs.PrintDefaults()
		fs.SetOutput(oldOutput)
		fmt.Fprintln(b)
		fmt.Fprintln(b, "Examples:")
		fmt.Fprintln(b, "  # Default: write run summary JSON to stdout")
		fmt.Fprintln(b, "  radigest -fasta ref.fa -enzymes EcoRI,MseI")
		fmt.Fprintln(b, "  # Pipe FASTA in and write GFF to stdout")
		fmt.Fprintln(b, "  zcat ref.fa.gz | radigest -fasta - -enzymes EcoRI,MseI -gff - | bgzip > frag.gff3.gz")
		fmt.Fprintln(b, "  # Single digest (EcoRI) to GFF file")
		fmt.Fprintln(b, "  radigest -fasta ref.fa -enzymes EcoRI -gff out.gff3")
		fmt.Fprintln(b, "  # Double digest with hard size selection + JSON summary file")
		fmt.Fprintln(b, "  radigest -fasta ref.fa -enzymes EcoRI,MseI -min 100 -max 800 -json run.json")
		fmt.Fprintln(b, "  # Double digest with soft-window scoring and broad per-fragment TSV")
		fmt.Fprintln(b, "  radigest -fasta ref.fa -enzymes PstI,MspI -min 250 -max 500 -score-min 1 -score-max 1000 -size-model soft-window -size-edge-sd 25 -fragments-tsv fragments.tsv -json run.json")
		fmt.Fprintln(b, "  # Simulate a 10 Mb genome at 42% GC and digest (chromosome name is always chr1)")
		fmt.Fprintln(b, "  radigest -sim-len 10000000 -sim-gc 0.42 -enzymes EcoRI,MseI -gff out.gff3")
		fmt.Fprint(stderr, b.String())
	}

	if err := fs.Parse(args); err != nil {
		if errors.Is(err, flag.ErrHelp) {
			return nil
		}
		return usageError{err: err}
	}

	gffOutputPath := normalizeOutputPath(*gffPath)
	fragmentsTSVOutputPath := normalizeOutputPath(*fragmentsTSVPath)
	fragmentsFASTAOutputPath := normalizeOutputPath(*fragmentsFASTAPath)
	jsonOutputPath := normalizeOutputPath(*jsonPath)
	if !anyFlagSet(fs, "gff", "fragments-tsv", "fragments-fasta", "json") {
		jsonOutputPath = "-"
	}

	if *showVer {
		fmt.Fprintf(stdout, "radigest %s\n", version)
		return nil
	}
	if *listEns {
		names := make([]string, 0, len(enzyme.DB))
		for name := range enzyme.DB {
			names = append(names, name)
		}
		sort.Strings(names)
		for _, n := range names {
			fmt.Fprintln(stdout, n)
		}
		return nil
	}

	// ---- validation ---------------------------------------------------------
	if *enzFlag == "" {
		fs.Usage()
		return usageError{err: errors.New("-enzymes is required")}
	}
	if (*fastaPath == "" && *simLen <= 0) || (*fastaPath != "" && *simLen > 0) {
		fs.Usage()
		return usageError{err: errors.New("use exactly one of -fasta or -sim-len")}
	}
	if err := validatePositiveThreads(*threads); err != nil {
		return err
	}
	if *minLen > *maxLen {
		return fmt.Errorf("invalid range: -min (%d) > -max (%d)", *minLen, *maxLen)
	}
	if *simLen > 0 {
		if err := validateSimGC(*simGC); err != nil {
			return err
		}
	}

	if err := validateOutputSelection(gffOutputPath, fragmentsTSVOutputPath, fragmentsFASTAOutputPath, jsonOutputPath); err != nil {
		return err
	}
	if err := validateOutputPaths(*fastaPath, gffOutputPath, fragmentsTSVOutputPath, fragmentsFASTAOutputPath, jsonOutputPath, *fastaPath != ""); err != nil {
		return err
	}

	scoreMin := *scoreMinFlag
	if scoreMin < 0 {
		scoreMin = *minLen
	}
	scoreMax := *scoreMaxFlag
	if scoreMax < 0 {
		scoreMax = *maxLen
	}
	selector, err := sizeselect.New(sizeselect.Config{
		Model:    sizeselect.Model(*sizeModel),
		Min:      *minLen,
		Max:      *maxLen,
		ScoreMin: scoreMin,
		ScoreMax: scoreMax,
		Mean:     *sizeMean,
		SD:       *sizeSD,
		EdgeSD:   *sizeEdgeSD,
	})
	if err != nil {
		return err
	}

	// Digest the union of the hard output window and the broader scoring window.
	// Optional writers decide which fragments are serialized to artifact outputs.
	digestMin := minInt(*minLen, scoreMin)
	digestMax := maxInt(*maxLen, scoreMax)

	// ---- compile enzymes ----------------------------------------------------
	ens, enzymeNames, err := parseEnzymes(*enzFlag)
	if err != nil {
		return err
	}
	plan, err := digest.TryNewPlanWithOptions(ens, digest.Options{
		AllowSame:   *allowSame,
		StrictCuts:  *strictCuts,
		IncludeEnds: *includeEnds,
	})
	if err != nil {
		return usageError{err: err}
	}

	// ---- start writers -------------------------------------------------------
	writer, err := collector.NewWriterTo(gffOutputPath, stdout)
	if err != nil {
		return fmt.Errorf("collector: %w", err)
	}
	fragWriter, err := fragmenttsv.NewTo(fragmentsTSVOutputPath, stdout)
	if err != nil {
		return fmt.Errorf("fragments tsv: %w", err)
	}
	fragFASTAWriter, err := fragmentfasta.NewTo(fragmentsFASTAOutputPath, stdout)
	if err != nil {
		return fmt.Errorf("fragments fasta: %w", err)
	}
	wantFragmentFASTA := fragmentsFASTAOutputPath != ""

	// ---- worker pool --------------------------------------------------------
	type job struct {
		idx int
		rec fasta.Record
	}
	jobs := make(chan job, *threads)
	results := make(chan digestResult, *threads)
	var wg sync.WaitGroup
	for i := 0; i < *threads; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for j := range jobs {
				fragCh := make(chan digest.Fragment, 64)
				errCh := make(chan error, 1)
				var seq []byte
				if wantFragmentFASTA {
					seq = j.rec.Seq
				}
				results <- digestResult{idx: j.idx, chr: j.rec.ID, seq: seq, frags: fragCh, errors: errCh}

				err := plan.DigestEach(j.rec.Seq, digestMin, digestMax, func(fr digest.Fragment) error {
					fragCh <- fr
					return nil
				})
				close(fragCh)
				errCh <- err
				close(errCh)
			}
		}()
	}
	go func() {
		wg.Wait()
		close(results)
	}()

	// ---- source: FASTA or synthetic ----------------------------------------
	resolvedSimSeed := int64(0)
	if *simLen > 0 {
		resolvedSimSeed = sim.ResolveSeed(*simSeed)
	}

	faCh := make(chan fasta.Record)
	sourceErrCh := make(chan error, 1)
	go func() {
		if *simLen > 0 {
			seq := sim.Make(*simLen, *simGC, resolvedSimSeed) // chr1
			faCh <- fasta.Record{ID: "chr1", Seq: seq}
			close(faCh)
			sourceErrCh <- nil
			return
		}
		sourceErrCh <- fasta.StreamFrom(*fastaPath, stdin, faCh)
	}()
	go func() {
		idx := 0
		for rec := range faCh {
			jobs <- job{idx: idx, rec: rec}
			idx++
		}
		close(jobs)
	}()

	// ---- wait + finalize ----------------------------------------------------
	sizeStats, streamErr := writeResultStreamsScoredTo(writer, fragWriter, fragFASTAWriter, selector, results, *verbose, stderr)
	stats, closeErr := writer.Close()
	fragCloseErr := fragWriter.Close()
	fragFASTACloseErr := fragFASTAWriter.Close()
	if streamErr != nil {
		return fmt.Errorf("digest/write: %w", streamErr)
	}
	sourceErr := <-sourceErrCh
	if sourceErr != nil {
		return fmt.Errorf("fasta stream: %w", sourceErr)
	}
	if closeErr != nil {
		return fmt.Errorf("collector: %w", closeErr)
	}
	if fragCloseErr != nil {
		return fmt.Errorf("fragments tsv: %w", fragCloseErr)
	}
	if fragFASTACloseErr != nil {
		return fmt.Errorf("fragments fasta: %w", fragFASTACloseErr)
	}

	fmt.Fprintf(stderr, "Fragments kept: %d\nBases covered: %d\nChromosomes: %d\n",
		stats.TotalFragments, stats.TotalBases, len(stats.PerChr))
	if jsonOutputPath != "" {
		summary := buildRunSummary(runSummaryInput{
			Args:               args,
			Enzymes:            enzymeNames,
			FastaPath:          *fastaPath,
			SimLen:             *simLen,
			SimGC:              *simGC,
			SimSeedRequested:   *simSeed,
			SimSeedResolved:    resolvedSimSeed,
			MinLen:             *minLen,
			MaxLen:             *maxLen,
			Threads:            *threads,
			AllowSame:          *allowSame,
			StrictCuts:         *strictCuts,
			IncludeEnds:        *includeEnds,
			SelectorConfig:     selector.Config(),
			JSONPath:           jsonOutputPath,
			GFFPath:            gffOutputPath,
			FragmentsTSVPath:   fragmentsTSVOutputPath,
			FragmentsFASTAPath: fragmentsFASTAOutputPath,
			SizeSelection:      sizeStats,
			Stats:              stats,
		})
		if err := writeSummaryJSONTo(jsonOutputPath, summary, stdout); err != nil {
			return fmt.Errorf("write json: %w", err)
		}
	}
	return nil
}

func writeResultStreams(w *collector.Writer, results <-chan digestResult, verbose bool) error {
	return writeResultStreamsTo(w, results, verbose, os.Stderr)
}

func writeResultStreamsTo(w *collector.Writer, results <-chan digestResult, verbose bool, stderr io.Writer) error {
	pending := make(map[int]digestResult)
	next := 0

	for results != nil || len(pending) > 0 {
		if r, ok := pending[next]; ok {
			cs, writeErr := w.WriteStream(r.chr, r.frags)
			digestErr := <-r.errors
			delete(pending, next)
			next++

			if verbose {
				fmt.Fprintf(stderr, "chr=%s fragments=%d\n", r.chr, cs.Fragments)
			}
			if writeErr != nil {
				return fmt.Errorf("write fragments for %s: %w", r.chr, writeErr)
			}
			if digestErr != nil {
				return fmt.Errorf("digest %s: %w", r.chr, digestErr)
			}
			continue
		}

		if results == nil {
			return fmt.Errorf("missing digest result for chromosome index %d", next)
		}
		r, ok := <-results
		if !ok {
			results = nil
			continue
		}
		pending[r.idx] = r
	}
	return nil
}

func writeResultStreamsScored(w *collector.Writer, tsv *fragmenttsv.Writer, fastaWriter *fragmentfasta.Writer, selector sizeselect.Selector, results <-chan digestResult, verbose bool) (sizeselect.Stats, error) {
	return writeResultStreamsScoredTo(w, tsv, fastaWriter, selector, results, verbose, os.Stderr)
}

func writeResultStreamsScoredTo(w *collector.Writer, tsv *fragmenttsv.Writer, fastaWriter *fragmentfasta.Writer, selector sizeselect.Selector, results <-chan digestResult, verbose bool, stderr io.Writer) (sizeselect.Stats, error) {
	pending := make(map[int]digestResult)
	next := 0
	stats := sizeselect.NewStats(selector)

	for results != nil || len(pending) > 0 {
		if r, ok := pending[next]; ok {
			cs, writeErr := writeScoredChromosome(w, tsv, fastaWriter, selector, &stats, r.chr, r.seq, r.frags)
			digestErr := <-r.errors
			delete(pending, next)
			next++

			if verbose {
				fmt.Fprintf(stderr, "chr=%s fragments=%d scored=%d\n", r.chr, cs.Fragments, stats.RawFragmentsScored)
			}
			if writeErr != nil {
				return stats, fmt.Errorf("write fragments for %s: %w", r.chr, writeErr)
			}
			if digestErr != nil {
				return stats, fmt.Errorf("digest %s: %w", r.chr, digestErr)
			}
			continue
		}

		if results == nil {
			return stats, fmt.Errorf("missing digest result for chromosome index %d", next)
		}
		r, ok := <-results
		if !ok {
			results = nil
			continue
		}
		pending[r.idx] = r
	}
	return stats, nil
}

func writeScoredChromosome(w *collector.Writer, tsv *fragmenttsv.Writer, fastaWriter *fragmentfasta.Writer, selector sizeselect.Selector, stats *sizeselect.Stats, chr string, seq []byte, frags <-chan digest.Fragment) (collector.ChrStats, error) {
	var local collector.ChrStats
	var firstErr error
	ordinal := 1

	for fr := range frags {
		length := fr.End - fr.Start
		hardKept := selector.InHardWindow(length)
		if hardKept {
			stats.AddHardKept(length)
		}
		if selector.InScoreRange(length) {
			weight := selector.Weight(length)
			stats.AddScored(length, weight)
			if firstErr == nil {
				if err := tsv.Write(chr, fr, hardKept, weight); err != nil {
					firstErr = err
				}
			}
		}
		if hardKept {
			if firstErr == nil {
				if err := w.WriteFragment(chr, ordinal, fr); err != nil {
					firstErr = err
				} else if err := fastaWriter.Write(chr, ordinal, fr, seq); err != nil {
					firstErr = err
				} else {
					local.Fragments++
					local.Bases += length
				}
			}
			ordinal++
		}
	}
	return local, firstErr
}

type runSummaryInput struct {
	Args               []string
	Enzymes            []string
	FastaPath          string
	SimLen             int
	SimGC              float64
	SimSeedRequested   int64
	SimSeedResolved    int64
	MinLen             int
	MaxLen             int
	Threads            int
	AllowSame          bool
	StrictCuts         bool
	IncludeEnds        bool
	SelectorConfig     sizeselect.Config
	JSONPath           string
	GFFPath            string
	FragmentsTSVPath   string
	FragmentsFASTAPath string
	SizeSelection      sizeselect.Stats
	Stats              collector.Stats
}

func buildRunSummary(in runSummaryInput) runSummary {
	params := parameterSummary{
		MinLength:   in.MinLen,
		MaxLength:   in.MaxLen,
		ScoreMin:    in.SelectorConfig.ScoreMin,
		ScoreMax:    in.SelectorConfig.ScoreMax,
		SizeModel:   string(in.SelectorConfig.Model),
		Threads:     in.Threads,
		AllowSame:   in.AllowSame,
		StrictCuts:  in.StrictCuts,
		IncludeEnds: in.IncludeEnds,
	}
	switch in.SelectorConfig.Model {
	case sizeselect.ModelNormal:
		params.SizeMean = in.SelectorConfig.Mean
		params.SizeSD = in.SelectorConfig.SD
	case sizeselect.ModelTriangular:
		params.SizeMean = in.SelectorConfig.Mean
	case sizeselect.ModelSoftWindow:
		params.SizeEdgeSD = in.SelectorConfig.EdgeSD
	}

	input := inputSummary{
		Source: "fasta",
		FASTA:  in.FastaPath,
	}
	warnings := []string{}
	if in.SimLen > 0 {
		simGC := in.SimGC
		simSeedRequested := in.SimSeedRequested
		simSeedResolved := in.SimSeedResolved
		input = inputSummary{
			Source:           "simulation",
			SimLength:        in.SimLen,
			SimGC:            &simGC,
			SimSeedRequested: &simSeedRequested,
			SimSeedResolved:  &simSeedResolved,
		}
		if in.SimSeedRequested == 0 {
			warnings = append(warnings, "-sim-seed 0 requested a time-based seed; resolved seed is recorded in input.sim_seed_resolved")
		}
	}
	if in.Stats.TotalFragments == 0 {
		warnings = append(warnings, "no fragments passed the hard size-selection window")
	}

	command := make([]string, 0, len(in.Args)+1)
	command = append(command, "radigest")
	command = append(command, in.Args...)

	outputs := outputSummary{
		JSON:           in.JSONPath,
		GFF:            in.GFFPath,
		FragmentsTSV:   in.FragmentsTSVPath,
		FragmentsFASTA: in.FragmentsFASTAPath,
	}

	return runSummary{
		SchemaVersion:   summarySchemaVersion,
		RadigestVersion: version,
		Command:         command,
		Input:           input,
		Parameters:      params,
		Outputs:         outputs,
		Warnings:        warnings,

		Enzymes:        in.Enzymes,
		MinLength:      in.MinLen,
		MaxLength:      in.MaxLen,
		GFF:            in.GFFPath,
		FragmentsTSV:   in.FragmentsTSVPath,
		FragmentsFASTA: in.FragmentsFASTAPath,
		SizeSelection:  in.SizeSelection,
		Stats:          in.Stats,
	}
}

func anyFlagSet(fs *flag.FlagSet, names ...string) bool {
	wanted := make(map[string]struct{}, len(names))
	for _, name := range names {
		wanted[name] = struct{}{}
	}
	found := false
	fs.Visit(func(f *flag.Flag) {
		if _, ok := wanted[f.Name]; ok {
			found = true
		}
	})
	return found
}

func normalizeOutputPath(path string) string {
	return strings.TrimSpace(path)
}

func writeSummaryJSON(path string, summary runSummary) error {
	return writeSummaryJSONTo(path, summary, os.Stdout)
}

func writeSummaryJSONTo(path string, summary runSummary, stdout io.Writer) error {
	if path == "" {
		return nil
	}
	if path == "-" {
		return json.NewEncoder(stdout).Encode(summary)
	}

	f, err := os.Create(path)
	if err != nil {
		return err
	}
	if err := json.NewEncoder(f).Encode(summary); err != nil {
		_ = f.Close()
		return err
	}
	return f.Close()
}

func minInt(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func maxInt(a, b int) int {
	if a > b {
		return a
	}
	return b
}
