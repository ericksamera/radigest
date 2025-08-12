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

func produceChunks(faCh <-chan fasta.Record, jobs chan<- fasta.Record, chunkSz int) {
	defer close(jobs)
	for rec := range faCh {
		seq := rec.Seq
		n := len(seq)
		for from := 0; from < n; from += chunkSz {
			to := from + chunkSz
			if to > n {
				to = n
			}
			jobs <- fasta.Record{
				ID:  rec.ID,
				Seq: seq[from:to],
			}
		}
	}
}

func main() {
	// ---- CLI flags ----------------------------------------------------------
	fastaPath := flag.String("fasta", "", "reference FASTA file (required)")
	enzFlag := flag.String("enzymes", "", "comma-separated enzyme names (≥2, first two form the AB pair)")
	minLen := flag.Int("min", 0, "minimum fragment length (bp)")
	maxLen := flag.Int("max", 1<<30, "maximum fragment length (bp)")
	gffPath := flag.String("gff", "fragments.gff3", "output GFF3 file")
	jsonPath := flag.String("json", "", "optional: write run summary JSON here")
	chunkSz := flag.Int("chunk", 8<<20, "chunk size (bp) sent to each worker")
	threads := flag.Int("threads", runtime.NumCPU(), "number of worker goroutines")
	showVer := flag.Bool("version", false, "print version and exit")
	listEns := flag.Bool("list-enzymes", false, "list available enzyme names and exit")

	flag.Usage = func() {
		b := &strings.Builder{}
		fmt.Fprintln(b, "radigest — in-silico double-digest and GFF3 fragment export")
		fmt.Fprintln(b)
		fmt.Fprintln(b, "Usage:")
		fmt.Fprintln(b, "  radigest -fasta <ref.fa> -enzymes <E1,E2[,E3...]> [options]")
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
		fmt.Fprintln(b, "  # Basic EcoRI/MseI digest to GFF3")
		fmt.Fprintln(b, "  radigest -fasta ref.fa -enzymes EcoRI,MseI -gff out.gff3")
		fmt.Fprintln(b, "  # Restrict fragment size and emit JSON summary")
		fmt.Fprintln(b, "  radigest -fasta ref.fa -enzymes EcoRI,MseI -min 100 -max 800 -json run.json")
		fmt.Fprintln(b, "  # See supported enzymes")
		fmt.Fprintln(b, "  radigest -list-enzymes")
		fmt.Fprintln(os.Stderr, b.String())
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
		fmt.Fprintln(os.Stderr, "error: flags -fasta and -enzymes are required\n")
		flag.Usage()
		os.Exit(2)
	}

	// ---- build enzyme slice -------------------------------------------------
	var ens []enzyme.Enzyme
	for _, name := range strings.Split(*enzFlag, ",") {
		e, ok := enzyme.DB[name]
		if !ok {
			log.Fatalf("unknown enzyme %q", name)
		}
		ens = append(ens, e)
	}
	if len(ens) < 2 {
		log.Fatal("need at least two enzymes")
	}

	// ---- start collector ----------------------------------------------------
	cIn, done, err := collector.New(*gffPath)
	if err != nil {
		log.Fatalf("collector: %v", err)
	}

	// ---- worker pool --------------------------------------------------------
	jobs := make(chan fasta.Record, *threads)
	var wg sync.WaitGroup
	for i := 0; i < *threads; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for rec := range jobs {
				frags := digest.Digest(rec.Seq, ens, *minLen, *maxLen)
				if len(frags) != 0 {
					cIn <- collector.Msg{Chr: rec.ID, Frags: frags}
				}
			}
		}()
	}

	// ---- stream FASTA into jobs --------------------------------------------
	faCh := make(chan fasta.Record)
	go func() {
		if err := fasta.Stream(*fastaPath, faCh); err != nil {
			log.Fatalf("fasta stream: %v", err)
		}
		// NOTE: assume fasta.Stream closes faCh when it returns.
	}()

	// single consumer / producer path
	go produceChunks(faCh, jobs, *chunkSz)

	// wait for workers, finish collector
	wg.Wait()  // jobs closed by produceChunks
	close(cIn) // tell collector to finish

	// ---- summary ------------------------------------------------------------
	stats := <-done
	fmt.Printf("Fragments kept: %d\nBases covered: %d\nChromosomes: %d\n",
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
