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

	"radigest/internal/collector"
	"radigest/internal/digest"
	"radigest/internal/enzyme"
	"radigest/internal/fasta"
	"radigest/internal/sim"
)

var (
	version = "v1.3.1"
)

func main() {
	// ---- CLI flags ----------------------------------------------------------
	fastaPath := flag.String("fasta", "", "reference FASTA file")
	enzFlag := flag.String("enzymes", "", "comma-separated enzyme names (≥1; first two form the AB pair if present)")
	minLen := flag.Int("min", 1, "minimum fragment length (bp)")
	maxLen := flag.Int("max", 1<<30, "maximum fragment length (bp)")
	gffPath := flag.String("gff", "fragments.gff3", "output GFF3 file (path or '-' for stdout)")
	jsonPath := flag.String("json", "", "optional: write run summary JSON here")
	threads := flag.Int("threads", runtime.NumCPU(), "number of worker goroutines")
	verbose := flag.Bool("v", false, "verbose progress to stderr")
	showVer := flag.Bool("version", false, "print version and exit")
	listEns := flag.Bool("list-enzymes", false, "list available enzyme names and exit")

	// double-digest behavior & validation
	allowSame := flag.Bool("allow-same", false, "double digest: also keep AA/BB neighbors (default AB/BA only)")
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
		fmt.Fprintln(b, "radigest — in-silico single/double digest and GFF3 fragment export")
		fmt.Fprintln(b)
		fmt.Fprintln(b, "Usage:")
		fmt.Fprintln(b, "  radigest -fasta <ref.fa|-> -enzymes <E1[,E2[,E3...]]> [options]")
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
		fmt.Fprintln(b, "  # Double digest with size selection + JSON summary")
		fmt.Fprintln(b, "  radigest -fasta ref.fa -enzymes EcoRI,MseI -min 100 -max 800 -json run.json")
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
	if *minLen > *maxLen {
		log.Fatalf("invalid range: -min (%d) > -max (%d)", *minLen, *maxLen)
	}

	// ---- compile enzymes ----------------------------------------------------
	var ens []enzyme.Enzyme
	for _, name := range strings.Split(*enzFlag, ",") {
		n := strings.TrimSpace(name)
		e, ok := enzyme.DB[n]
		if !ok {
			log.Fatalf("unknown enzyme %q", name)
		}
		ens = append(ens, e)
	}
	if len(ens) >= 2 && ens[0].Name == ens[1].Name {
		log.Fatalf("first two enzymes must differ (got %s,%s)", ens[0].Name, ens[1].Name)
	}
	plan := digest.NewPlanWithOptions(ens, digest.Options{
		AllowSame:  *allowSame,
		StrictCuts: *strictCuts,
	})

	// ---- start collector ----------------------------------------------------
	cIn, done, err := collector.New(*gffPath)
	if err != nil {
		log.Fatalf("collector: %v", err)
	}

	// ---- worker pool --------------------------------------------------------
	type job struct {
		idx int
		rec fasta.Record
	}
	jobs := make(chan job, *threads)
	var wg sync.WaitGroup
	for i := 0; i < *threads; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for j := range jobs {
				frags := plan.Digest(j.rec.Seq, *minLen, *maxLen)
				cIn <- collector.Msg{Idx: j.idx, Chr: j.rec.ID, Frags: frags}
				if *verbose {
					fmt.Fprintf(os.Stderr, "chr=%s fragments=%d\n", j.rec.ID, len(frags))
				}
			}
		}()
	}

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
	wg.Wait()
	close(cIn)

	stats := <-done
	fmt.Fprintf(os.Stderr, "Fragments kept: %d\nBases covered: %d\nChromosomes: %d\n",
		stats.TotalFragments, stats.TotalBases, len(stats.PerChr))
	if *jsonPath != "" {
		out := struct {
			Enzymes   []string `json:"enzymes"`
			MinLength int      `json:"min_length"`
			MaxLength int      `json:"max_length"`
			collector.Stats
		}{
			Enzymes:   strings.Split(*enzFlag, ","),
			MinLength: *minLen,
			MaxLength: *maxLen,
			Stats:     stats,
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
