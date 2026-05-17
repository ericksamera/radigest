package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"log"
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

type digestResult struct {
	idx    int
	chr    string
	seq    []byte
	frags  <-chan digest.Fragment
	errors <-chan error
}

func main() {
	// ---- CLI flags ----------------------------------------------------------
	fastaPath := flag.String("fasta", "", "reference FASTA file")
	enzFlag := flag.String("enzymes", "", "comma-separated enzyme names (one or two; two form the AB pair)")
	minLen := flag.Int("min", 1, "minimum fragment length (bp) for hard-selected GFF output")
	maxLen := flag.Int("max", 1<<30, "maximum fragment length (bp) for hard-selected GFF output")
	gffPath := flag.String("gff", "fragments.gff3", "output GFF3 file for hard-selected fragments (path or '-' for stdout)")
	fragmentsTSVPath := flag.String("fragments-tsv", "fragments.tsv", "per-fragment TSV for score-range fragments; empty string disables")
	fragmentsFASTAPath := flag.String("fragments-fasta", "", "FASTA output for hard-selected fragments; empty string disables")
	jsonPath := flag.String("json", "", "optional: write run summary JSON here")
	threads := flag.Int("threads", runtime.NumCPU(), "number of worker goroutines")
	verbose := flag.Bool("v", false, "verbose progress to stderr")
	showVer := flag.Bool("version", false, "print version and exit")
	listEns := flag.Bool("list-enzymes", false, "list available enzyme names and exit")

	// size-selection scoring
	scoreMinFlag := flag.Int("score-min", -1, "minimum fragment length included in fragments TSV and size-selection stats; default -min")
	scoreMaxFlag := flag.Int("score-max", -1, "maximum fragment length included in fragments TSV and size-selection stats; default -max")
	sizeModel := flag.String("size-model", "hard", "size-selection model: hard, normal, triangular, or soft-window")
	sizeMean := flag.Float64("size-mean", 0, "target/peak insert length for normal/triangular models; default midpoint of -min/-max")
	sizeSD := flag.Float64("size-sd", 35, "standard deviation for -size-model normal")
	sizeEdgeSD := flag.Float64("size-edge-sd", 25, "edge softness for -size-model soft-window")

	// digest behavior & validation
	allowSame := flag.Bool("allow-same", false, "double digest: also keep AA/BB neighbors (default AB/BA only)")
	includeEnds := flag.Bool("include-ends", false, "also emit terminal fragments from chromosome/contig ends to the nearest cut")
	strictCuts := flag.Bool("strict-cuts", false, "error if an enzyme lacks a caret and CutIndex==0 (no mid-site fallback)")

	// synthetic genome flags
	simLen := flag.Int("sim-len", 0, "synthesize a single-chromosome genome of this length (bp) instead of reading -fasta")
	simGC := flag.Float64("sim-gc", 0.50, "target GC fraction in [0,1] for -sim-len")
	simSeed := flag.Int64("sim-seed", 1, "PRNG seed for -sim-len (0 ⇒ time-based)")

	flag.Usage = func() {
		b := &strings.Builder{}
		fmt.Fprintln(b, "Author:  Erick Samera (erick.samera@kpu.ca)")
		fmt.Fprintln(b, "License: MIT")
		fmt.Fprintln(b, "Version:", version)
		fmt.Fprintln(b)
		fmt.Fprintln(b, "radigest — in-silico single/double digest and GFF3/TSV fragment export")
		fmt.Fprintln(b)
		fmt.Fprintln(b, "Usage:")
		fmt.Fprintln(b, "  radigest -fasta <ref.fa|-> -enzymes <E1[,E2]> [options]")
		fmt.Fprintln(b, "  radigest -sim-len <bp> -sim-gc <0..1> -enzymes <E1[,E2]> [options]")
		fmt.Fprintln(b)
		fmt.Fprintln(b, "Required flags:")
		fmt.Fprintln(b, "  -enzymes, and exactly one of -fasta or -sim-len")
		fmt.Fprintln(b)
		fmt.Fprintln(b, "Options:")
		flag.CommandLine.SetOutput(b)
		flag.PrintDefaults()
		flag.CommandLine.SetOutput(os.Stderr)
		fmt.Fprintln(b)
		fmt.Fprintln(b, "Examples:")
		fmt.Fprintln(b, "  # Pipe FASTA in and write GFF to stdout")
		fmt.Fprintln(b, "  zcat ref.fa.gz | radigest -fasta - -enzymes EcoRI,MseI -gff - | bgzip > frag.gff3.gz")
		fmt.Fprintln(b, "  # Single digest (EcoRI) to file")
		fmt.Fprintln(b, "  radigest -fasta ref.fa -enzymes EcoRI -gff out.gff3")
		fmt.Fprintln(b, "  # Double digest with hard size selection + JSON summary")
		fmt.Fprintln(b, "  radigest -fasta ref.fa -enzymes EcoRI,MseI -min 100 -max 800 -json run.json")
		fmt.Fprintln(b, "  # Double digest with soft-window scoring and broad per-fragment TSV")
		fmt.Fprintln(b, "  radigest -fasta ref.fa -enzymes PstI,MspI -min 250 -max 500 -score-min 1 -score-max 1000 -size-model soft-window -size-edge-sd 25 -fragments-tsv fragments.tsv -json run.json")
		fmt.Fprintln(b, "  # Simulate a 10 Mb genome at 42% GC and digest (chromosome name is always chr1)")
		fmt.Fprintln(b, "  radigest -sim-len 10000000 -sim-gc 0.42 -enzymes EcoRI,MseI -gff out.gff3")
		fmt.Fprint(os.Stderr, b.String())
	}

	flag.Parse()

	if *showVer {
		fmt.Printf("radigest %s\n", version)
		return
	}
	if *listEns {
		names := make([]string, 0, len(enzyme.DB))
		for name := range enzyme.DB {
			names = append(names, name)
		}
		sort.Strings(names)
		for _, n := range names {
			fmt.Println(n)
		}
		return
	}

	// ---- validation ---------------------------------------------------------
	if *enzFlag == "" {
		fmt.Fprintln(os.Stderr, "error: -enzymes is required")
		flag.Usage()
		os.Exit(2)
	}
	if (*fastaPath == "" && *simLen <= 0) || (*fastaPath != "" && *simLen > 0) {
		fmt.Fprintln(os.Stderr, "error: use exactly one of -fasta or -sim-len")
		flag.Usage()
		os.Exit(2)
	}
	if err := validatePositiveThreads(*threads); err != nil {
		log.Fatal(err)
	}
	if *minLen > *maxLen {
		log.Fatalf("invalid range: -min (%d) > -max (%d)", *minLen, *maxLen)
	}
	if *simLen > 0 {
		if err := validateSimGC(*simGC); err != nil {
			log.Fatal(err)
		}
	}

	if err := validateOutputPaths(*fastaPath, *gffPath, *fragmentsTSVPath, *fragmentsFASTAPath, *jsonPath, *fastaPath != ""); err != nil {
		log.Fatal(err)
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
		log.Fatal(err)
	}

	// Digest the union of the hard GFF window and the broader scoring window.
	// The writer decides which fragments go to GFF and which go to TSV/stats.
	digestMin := minInt(*minLen, scoreMin)
	digestMax := maxInt(*maxLen, scoreMax)

	// ---- compile enzymes ----------------------------------------------------
	ens, enzymeNames, err := parseEnzymes(*enzFlag)
	if err != nil {
		log.Fatal(err)
	}
	plan := digest.NewPlanWithOptions(ens, digest.Options{
		AllowSame:   *allowSame,
		StrictCuts:  *strictCuts,
		IncludeEnds: *includeEnds,
	})

	// ---- start writers -------------------------------------------------------
	writer, err := collector.NewWriter(*gffPath)
	if err != nil {
		log.Fatalf("collector: %v", err)
	}
	fragWriter, err := fragmenttsv.New(*fragmentsTSVPath)
	if err != nil {
		log.Fatalf("fragments tsv: %v", err)
	}
	fragFASTAWriter, err := fragmentfasta.New(*fragmentsFASTAPath)
	if err != nil {
		log.Fatalf("fragments fasta: %v", err)
	}
	wantFragmentFASTA := strings.TrimSpace(*fragmentsFASTAPath) != ""

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
	faCh := make(chan fasta.Record)
	go func() {
		if *simLen > 0 {
			seq := sim.Make(*simLen, *simGC, *simSeed) // chr1
			faCh <- fasta.Record{ID: "chr1", Seq: seq}
			close(faCh)
			return
		}
		if err := fasta.Stream(*fastaPath, faCh); err != nil {
			log.Fatalf("fasta stream: %v", err)
		}
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
	sizeStats, streamErr := writeResultStreamsScored(writer, fragWriter, fragFASTAWriter, selector, results, *verbose)
	stats, closeErr := writer.Close()
	fragCloseErr := fragWriter.Close()
	fragFASTACloseErr := fragFASTAWriter.Close()
	if streamErr != nil {
		log.Fatalf("digest/write: %v", streamErr)
	}
	if closeErr != nil {
		log.Fatalf("collector: %v", closeErr)
	}
	if fragCloseErr != nil {
		log.Fatalf("fragments tsv: %v", fragCloseErr)
	}
	if fragFASTACloseErr != nil {
		log.Fatalf("fragments fasta: %v", fragFASTACloseErr)
	}

	fmt.Fprintf(os.Stderr, "Fragments kept: %d\nBases covered: %d\nChromosomes: %d\n",
		stats.TotalFragments, stats.TotalBases, len(stats.PerChr))
	if *jsonPath != "" {
		out := struct {
			Enzymes        []string         `json:"enzymes"`
			MinLength      int              `json:"min_length"`
			MaxLength      int              `json:"max_length"`
			FragmentsTSV   string           `json:"fragments_tsv,omitempty"`
			FragmentsFASTA string           `json:"fragments_fasta,omitempty"`
			SizeSelection  sizeselect.Stats `json:"size_selection"`
			collector.Stats
		}{
			Enzymes:        enzymeNames,
			MinLength:      *minLen,
			MaxLength:      *maxLen,
			FragmentsTSV:   *fragmentsTSVPath,
			FragmentsFASTA: *fragmentsFASTAPath,
			SizeSelection:  sizeStats,
			Stats:          stats,
		}
		f, err := os.Create(*jsonPath)
		if err != nil {
			log.Fatalf("write json: %v", err)
		}
		if err := json.NewEncoder(f).Encode(out); err != nil {
			log.Fatalf("encode json: %v", err)
		}
		_ = f.Close()
	}
}

func writeResultStreams(w *collector.Writer, results <-chan digestResult, verbose bool) error {
	pending := make(map[int]digestResult)
	next := 0

	for results != nil || len(pending) > 0 {
		if r, ok := pending[next]; ok {
			cs, writeErr := w.WriteStream(r.chr, r.frags)
			digestErr := <-r.errors
			delete(pending, next)
			next++

			if verbose {
				fmt.Fprintf(os.Stderr, "chr=%s fragments=%d\n", r.chr, cs.Fragments)
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
				fmt.Fprintf(os.Stderr, "chr=%s fragments=%d scored=%d\n", r.chr, cs.Fragments, stats.RawFragmentsScored)
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
