package collector

import (
	"bufio"
	"fmt"
	"os"

	"radigest/internal/digest"
)

// Msg delivers a batch of fragments for one chromosome.
type Msg struct {
	Idx   int
	Chr   string
	Frags []digest.Fragment
}

type ChrStats struct {
	Fragments int `json:"fragments"`
	Bases     int `json:"bases"`
}

type Stats struct {
	TotalFragments int                 `json:"total_fragments"`
	TotalBases     int                 `json:"total_bases"`
	PerChr         map[string]ChrStats `json:"per_chromosome"`
}

// New starts the collector goroutine.
//   • send Msg values on the returned chan
//   • close the chan when workers are done
//   • read the final Stats from the second chan
func New(gffPath string) (chan<- Msg, <-chan Stats, error) {
	var f *os.File
	var err error
	if gffPath == "-" {
		f = os.Stdout
	} else {
		f, err = os.Create(gffPath)
		if err != nil {
			return nil, nil, err
		}
	}
	bw := bufio.NewWriter(f)

	// header is written once
	if _, err := bw.WriteString("##gff-version 3\n"); err != nil {
		f.Close()
		return nil, nil, err
	}

	in := make(chan Msg)
	out := make(chan Stats, 1)

	go func() {
		defer f.Close()
		defer close(out)

		stats := Stats{PerChr: make(map[string]ChrStats)}

		next := 0
		pending := make(map[int]Msg)

		write := func(msg Msg) {
			for i, fr := range msg.Frags {
				start := fr.Start + 1 // 1-based closed for GFF
				end := fr.End
				ln := end - fr.Start
				fmt.Fprintf(bw,
					"%s\tradigest\tfragment\t%d\t%d\t.\t+\t.\tID=%s_%d;Length=%d\n",
					msg.Chr, start, end, msg.Chr, i+1, ln)
				stats.TotalFragments++
				stats.TotalBases += ln
				cs := stats.PerChr[msg.Chr]
				cs.Fragments++
				cs.Bases += ln
				stats.PerChr[msg.Chr] = cs
			}
		}

		for msg := range in {
			pending[msg.Idx] = msg
			for {
				if m, ok := pending[next]; ok {
					write(m)
					delete(pending, next)
					next++
				} else {
					break
				}
			}
		}
		// defensive drain (should be empty)
		for ; len(pending) > 0; next++ {
			if m, ok := pending[next]; ok {
				write(m)
				delete(pending, next)
			} else {
				break
			}
		}

		if err := bw.Flush(); err != nil {
			fmt.Fprintf(os.Stderr, "collector flush: %v\n", err)
		}
		out <- stats
	}()

	return in, out, nil
}
