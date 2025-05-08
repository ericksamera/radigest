package main

import (
	"flag"
	"fmt"
	"log"
	"runtime"
	"strings"
	"sync"
	"encoding/json"
	"os"

	"radigest/internal/collector"
	"radigest/internal/digest"
	"radigest/internal/enzyme"
	"radigest/internal/fasta"
)

func main() {
	// ---- CLI flags ----------------------------------------------------------
	fastaPath := flag.String("fasta", "", "reference FASTA file (required)")
	enzFlag   := flag.String("enzymes", "", "comma-separated enzyme names (≥2, first two form the AB pair)")
	minLen    := flag.Int("min", 0, "minimum fragment length")
	maxLen    := flag.Int("max", 1<<30, "maximum fragment length")
	gffPath   := flag.String("gff",  "fragments.gff3", "output GFF3")
	jsonPath  := flag.String("json", "", "write run summary JSON")
	threads   := flag.Int("threads", runtime.NumCPU(), "worker goroutines")
	flag.Parse()

	if *fastaPath == "" || *enzFlag == "" {
		log.Fatal("flags --fasta and --enzymes are required")
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
	faCh := make(chan fasta.Record, 2)
	go func() {
		if err := fasta.Stream(*fastaPath, faCh); err != nil {
			log.Fatalf("fasta stream: %v", err)
		}
	}()

	for rec := range faCh {
		jobs <- rec
	}
	close(jobs) // no more work
	wg.Wait()   // workers done
	close(cIn)  // tell collector to finish

	// ---- summary ------------------------------------------------------------
	stats := <-done
	fmt.Printf("Fragments kept: %d\nBases covered: %d\nChromosomes: %d\n",
	    stats.TotalFragments, stats.TotalBases, len(stats.PerChr))	
	if *jsonPath != "" {
	    out := struct {
	        Enzymes     []string               `json:"enzymes"`
	        MinLength   int                    `json:"min_length"`
	        MaxLength   int                    `json:"max_length"`
	        collector.Stats
	    }{
	        Enzymes:   strings.Split(*enzFlag, ","),
	        MinLength: *minLen,
	        MaxLength: *maxLen,
	        Stats:     stats,
	    }
	    f, err := os.Create(*jsonPath)
	    if err != nil { log.Fatalf("write json: %v", err) }
	    if err := json.NewEncoder(f).Encode(out); err != nil {
	        log.Fatalf("encode json: %v", err)
	    }
	    f.Close()
	}
}