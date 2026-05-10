package collector

import (
	"bufio"
	"fmt"
	"os"

	"github.com/KPU-AGC/radigest/internal/digest"
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

// Writer serializes fragments to GFF3 and accumulates run statistics.
type Writer struct {
	f         *os.File
	bw        *bufio.Writer
	closeFile bool
	stats     Stats
}

// NewWriter opens gffPath, writes the GFF3 header, and returns a streaming
// writer. Use Close to flush and, when appropriate, close the underlying file.
func NewWriter(gffPath string) (*Writer, error) {
	var f *os.File
	closeFile := false
	if gffPath == "-" {
		f = os.Stdout
	} else {
		var err error
		f, err = os.Create(gffPath)
		if err != nil {
			return nil, err
		}
		closeFile = true
	}

	w := &Writer{
		f:         f,
		bw:        bufio.NewWriter(f),
		closeFile: closeFile,
		stats:     Stats{PerChr: make(map[string]ChrStats)},
	}
	if _, err := w.bw.WriteString("##gff-version 3\n"); err != nil {
		if closeFile {
			_ = f.Close()
		}
		return nil, err
	}
	return w, nil
}

func (w *Writer) WriteFragment(chr string, ordinal int, fr digest.Fragment) error {
	start := fr.Start + 1 // 1-based closed for GFF
	end := fr.End
	ln := end - fr.Start
	if _, err := fmt.Fprintf(w.bw,
		"%s\tradigest\tfragment\t%d\t%d\t.\t+\t.\tID=%s_%d;Length=%d\n",
		chr, start, end, chr, ordinal, ln); err != nil {
		return err
	}
	w.stats.TotalFragments++
	w.stats.TotalBases += ln
	cs := w.stats.PerChr[chr]
	cs.Fragments++
	cs.Bases += ln
	w.stats.PerChr[chr] = cs
	return nil
}

// WriteFragments writes one chromosome worth of fragments from a slice. It is
// retained for compatibility with the batch collector API.
func (w *Writer) WriteFragments(chr string, frags []digest.Fragment) error {
	for i, fr := range frags {
		if err := w.WriteFragment(chr, i+1, fr); err != nil {
			return err
		}
	}
	return nil
}

// WriteStream writes one chromosome worth of fragments from a channel. It drains
// the channel even after the first write error so upstream digest goroutines are
// not left blocked on sends.
func (w *Writer) WriteStream(chr string, frags <-chan digest.Fragment) (ChrStats, error) {
	var local ChrStats
	var firstErr error
	ordinal := 1
	for fr := range frags {
		if firstErr == nil {
			if err := w.WriteFragment(chr, ordinal, fr); err != nil {
				firstErr = err
			} else {
				local.Fragments++
				local.Bases += fr.End - fr.Start
			}
		}
		ordinal++
	}
	return local, firstErr
}

// Close flushes pending output, closes owned files, and returns accumulated
// statistics. Stdout is flushed but not closed.
func (w *Writer) Close() (Stats, error) {
	err := w.bw.Flush()
	if w.closeFile {
		if closeErr := w.f.Close(); err == nil {
			err = closeErr
		}
	}
	return w.stats, err
}

// New starts the collector goroutine.
//   - send Msg values on the returned chan
//   - close the chan when workers are done
//   - read the final Stats from the second chan
func New(gffPath string) (chan<- Msg, <-chan Stats, error) {
	w, err := NewWriter(gffPath)
	if err != nil {
		return nil, nil, err
	}

	in := make(chan Msg)
	out := make(chan Stats, 1)

	go func() {
		defer close(out)

		next := 0
		pending := make(map[int]Msg)

		write := func(msg Msg) {
			if err := w.WriteFragments(msg.Chr, msg.Frags); err != nil {
				fmt.Fprintf(os.Stderr, "collector write: %v\n", err)
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

		stats, err := w.Close()
		if err != nil {
			fmt.Fprintf(os.Stderr, "collector flush: %v\n", err)
		}
		out <- stats
	}()

	return in, out, nil
}
