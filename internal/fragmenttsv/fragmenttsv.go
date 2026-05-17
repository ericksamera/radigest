package fragmenttsv

import (
	"bufio"
	"fmt"
	"io"
	"os"

	"github.com/ericksamera/radigest/internal/digest"
)

// Writer emits per-fragment TSV rows for downstream modeling. A Writer created
// with an empty path is a no-op, which lets callers keep TSV output disabled
// without nil checks.
type Writer struct {
	bw       *bufio.Writer
	close    func() error
	disabled bool
}

// New opens path and writes the TSV header. Use an empty path to disable TSV
// output. Use "-" to write to stdout.
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
	w := &Writer{bw: bufio.NewWriter(sink), close: close}
	if _, err := w.bw.WriteString("chrom\tstart0\tend0\tlength\thard_kept\tsize_weight\n"); err != nil {
		if close != nil {
			_ = close()
		}
		return nil, err
	}
	return w, nil
}

// Write emits one scored fragment row. Coordinates are 0-based half-open.
func (w *Writer) Write(chr string, fr digest.Fragment, hardKept bool, sizeWeight float64) error {
	if w == nil || w.disabled {
		return nil
	}
	length := fr.End - fr.Start
	_, err := fmt.Fprintf(w.bw, "%s\t%d\t%d\t%d\t%t\t%.8g\n", chr, fr.Start, fr.End, length, hardKept, sizeWeight)
	return err
}

// Close flushes pending TSV output and closes owned files. Stdout is flushed but
// not closed. Disabled writers are no-ops.
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
