package collector

import (
	"bufio"
	"fmt"
	"os"

	"radigest/internal/digest"
)

// Msg delivers a batch of fragments for one chromosome.
type Msg struct {
	Chr   string
	Frags []digest.Fragment
}

// Stats is emitted after the input channel closes.
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
	f, err := os.Create(gffPath)
	if err != nil {
		return nil, nil, err
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

		for msg := range in {
			for i, fr := range msg.Frags {
				start := fr.Start + 1 // 1-based closed for GFF
				end := fr.End
				fmt.Fprintf(
					bw,
					"%s\tradigest\tfragment\t%d\t%d\t.\t+\t.\tID=%s_%d\n",
					msg.Chr, start, end, msg.Chr, i+1,
				)
				ln := end - fr.Start
				stats.TotalFragments++
				stats.TotalBases += ln
				cs := stats.PerChr[msg.Chr]
				cs.Fragments++
				cs.Bases += ln
				stats.PerChr[msg.Chr] = cs
			}
		}
		bw.Flush()
		out <- stats
	}()

	return in, out, nil
}
