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
)

var (
	version = "dev"
	commit  = "none"
	date    = "unknown"
)

func main() {
	// ---- CLI flags ----------------------------------------------------------
	fastaPath := flag.String("fasta", "", "reference FASTA file (required)")
	enzFlag := flag.String("enzymes", "", "comma-separated enzyme names (≥1; first two form the AB pair if present)")
	minLen := flag.Int("min", 0, "minimum fragment length (bp)")
	maxLen := flag.Int("max", 1<<30, "maximum fragment length (bp)")
	gffPath := flag.String("gff", "fragments.gff3", "output GFF3 file (path or '-' for stdout)")
	jsonPath := flag.String("json", "", "optional: write run summary JSON here")
	threads := flag.Int("threads", runtime.NumCPU(), "number of worker goroutines")
	verbose := flag.Bool("v", false, "verbose progress to stderr")
	showVer := flag.Bool("version", false, "print version and exit")
	listEns := flag.Bool("list-enzymes", false, "list available enzyme names and exit")

	flag.Usage = func() {
		b := &strings.Builder{}
		fmt.Fprintln(b, "radigest — in-silico single/double digest and GFF3 fragment export")
		fmt.Fprintln(b)
		fmt.Fprintln(b, "Usage:")
		fmt.Fprintln(b, "  radigest -fasta <ref.fa|-> -enzymes <E1[,E2[,E3...]]> [options]")
		fmt.Fprintln(b)
		fmt.Fprintln(b, "Required flags:")
		fmt.Fprintln(b, "  -fasta, -enzymes")
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
		fmt.Fprint(os.Stderr, b.String()) // avoid extra blank line
	}

	flag.Parse()

	if *showVer {
		fmt.Printf("radigest %s (commit %s, %s)\n", version, commit, date)
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
	if *fastaPath == "" || *enzFlag == "" {
		fmt.Fprintln(os.Stderr, "error: flags -fasta and -enzymes are required")
		flag.Usage()
		os.Exit(2)
	}

	// ---- build enzyme slice -------------------------------------------------
	var ens []enzyme.Enzyme
	for _, name := range strings.Split(*enzFlag, ",") {
		n := strings.TrimSpace(name) // forgive spaces
		e, ok := enzyme.DB[n]
		if !ok {
			log.Fatalf("unknown enzyme %q", name)
		}
		ens = append(ens, e)
	}
	if len(ens) < 1 {
		log.Fatal("need at least one enzyme")
	}
	if len(ens) >= 2 && ens[0].Name == ens[1].Name {
		log.Fatalf("first two enzymes must differ (got %s,%s)", ens[0].Name, ens[1].Name)
	}
	if *minLen > *maxLen {
		log.Fatalf("invalid range: -min (%d) > -max (%d)", *minLen, *maxLen)
	}

	// ---- compile plan once (perf) ------------------------------------------
	plan := digest.NewPlan(ens)

	// ---- start collector ----------------------------------------------------
	cIn, done, err := collector.New(*gffPath)
	if err != nil {
		log.Fatalf("collector: %v", err)
	}

	// ---- worker pool (deterministic order via job idx) ---------------------
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
				// send even if empty to advance deterministic ordering
				cIn <- collector.Msg{Idx: j.idx, Chr: j.rec.ID, Frags: frags}
				if *verbose {
					fmt.Fprintf(os.Stderr, "chr=%s fragments=%d\n", j.rec.ID, len(frags))
				}
			}
		}()
	}

	// ---- stream FASTA and feed jobs ----------------------------------------
	faCh := make(chan fasta.Record)
	go func() {
		if err := fasta.Stream(*fastaPath, faCh); err != nil {
			log.Fatalf("fasta stream: %v", err)
		}
		// Stream closes faCh when done.
	}()
	go func() {
		idx := 0
		for rec := range faCh {
			jobs <- job{idx: idx, rec: rec} // whole contig per job
			idx++
		}
		close(jobs)
	}()

	// wait for workers, finish collector
	wg.Wait()
	close(cIn)

	// ---- summary (to stderr so stdout stays pure GFF) ----------------------
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
