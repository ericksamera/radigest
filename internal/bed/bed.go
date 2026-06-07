// Package bed writes hard-kept digest fragments as BED6 records.
package bed

import (
	"bufio"
	"fmt"
	"io"
	"os"

	"github.com/ericksamera/radigest/internal/digest"
	"github.com/ericksamera/radigest/internal/gff"
)

// Writer emits BED6 rows for hard-kept fragments. A Writer created with an
// empty path is a no-op, which lets callers keep BED output disabled without
// nil checks.
type Writer struct {
	bw       *bufio.Writer
	close    func() error
	disabled bool
}

// New opens path for BED output. Use an empty path to disable BED output. Use
// "-" to write to stdout.
func New(path string) (*Writer, error) {
	return NewTo(path, os.Stdout)
}

// NewTo is like New, but writes "-" to stdout instead of os.Stdout.
func NewTo(path string, stdout io.Writer) (*Writer, error) {
	if path == "" {
		return &Writer{disabled: true}, nil
	}

	var sink io.Writer
	var close func() error
	if path == "-" {
		if stdout == nil {
			return nil, fmt.Errorf("stdout writer is nil")
		}
		sink = stdout
	} else {
		f, err := os.Create(path)
		if err != nil {
			return nil, err
		}
		sink = f
		close = f.Close
	}
	return &Writer{bw: bufio.NewWriter(sink), close: close}, nil
}

// Write emits one BED6 record. Coordinates are 0-based half-open, matching BED
// convention and the fragment TSV/FASTA metadata. The ordinal should match the
// corresponding saved fragment ordinal for the chromosome.
func (w *Writer) Write(chr string, ordinal int, fr digest.Fragment) error {
	if w == nil || w.disabled {
		return nil
	}
	if fr.Start < 0 || fr.End < fr.Start {
		return fmt.Errorf("BED: invalid fragment for %s_%d: start=%d end=%d", chr, ordinal, fr.Start, fr.End)
	}
	_, err := fmt.Fprintf(
		w.bw,
		"%s\t%d\t%d\t%s\t0\t+\n",
		bedChrom(chr),
		fr.Start,
		fr.End,
		fragmentID(chr, ordinal),
	)
	return err
}

// Close flushes pending BED output and closes owned files. Stdout is flushed
// but not closed. Disabled writers are no-ops.
func (w *Writer) Close() error {
	if w == nil || w.disabled {
		return nil
	}
	err := w.bw.Flush()
	if w.close != nil {
		if closeErr := w.close(); err == nil {
			err = closeErr
		}
	}
	return err
}

func bedChrom(chr string) string {
	if chr == "" {
		return "."
	}
	return chr
}

func fragmentID(chr string, ordinal int) string {
	if chr == "" {
		return fmt.Sprintf("frag%d", ordinal)
	}
	return fmt.Sprintf("%s_%d", gff.EscapeAttributeValue(chr), ordinal)
}
