package collector

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/ericksamera/radigest/internal/digest"
	"github.com/ericksamera/radigest/internal/gff"
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

// Writer optionally serializes fragments to GFF3 and accumulates run statistics.
// A Writer created with an empty path is stats-only and writes no GFF output.
type Writer struct {
	bw       *bufio.Writer
	close    func() error
	disabled bool
	stats    Stats
}

// NewWriter opens gffPath, writes the GFF3 header, and returns a streaming
// writer. Use an empty path for stats-only mode with GFF output disabled. Use
// Close to flush and, when appropriate, close the underlying file.
func NewWriter(gffPath string) (*Writer, error) {
	return NewWriterTo(gffPath, os.Stdout)
}

// NewWriterTo is like NewWriter, but writes "-" to stdout instead of os.Stdout.
func NewWriterTo(gffPath string, stdout io.Writer) (*Writer, error) {
	if strings.TrimSpace(gffPath) == "" {
		return &Writer{disabled: true, stats: Stats{PerChr: make(map[string]ChrStats)}}, nil
	}

	var sink io.Writer
	var close func() error
	if gffPath == "-" {
		if stdout == nil {
			return nil, fmt.Errorf("stdout writer is nil")
		}
		sink = stdout
	} else {
		f, err := os.Create(gffPath)
		if err != nil {
			return nil, err
		}
		sink = f
		close = f.Close
	}

	w := &Writer{
		bw:    bufio.NewWriter(sink),
		close: close,
		stats: Stats{PerChr: make(map[string]ChrStats)},
	}
	if _, err := w.bw.WriteString("##gff-version 3\n"); err != nil {
		if close != nil {
			_ = close()
		}
		return nil, err
	}
	return w, nil
}

// WriteFragment records one hard-kept fragment and, when GFF output is enabled,
// writes one GFF3 feature using the caller-provided per-chromosome ordinal.
func (w *Writer) WriteFragment(chr string, ordinal int, fr digest.Fragment) error {
	if w == nil {
		return nil
	}
	start := fr.Start + 1 // 1-based closed for GFF
	end := fr.End
	ln := end - fr.Start
	if !w.disabled {
		if _, err := fmt.Fprintf(w.bw,
			"%s\tradigest\tfragment\t%d\t%d\t.\t+\t.\t%s\n",
			gff.EscapeSeqID(chr), start, end, gff.FragmentAttributes(chr, ordinal, ln)); err != nil {
			return err
		}
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
// statistics. Stdout is flushed but not closed. Disabled writers are no-ops
// except for returning accumulated statistics.
func (w *Writer) Close() (Stats, error) {
	if w == nil {
		return Stats{}, nil
	}
	if w.disabled {
		return w.stats, nil
	}
	err := w.bw.Flush()
	if w.close != nil {
		if closeErr := w.close(); err == nil {
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
