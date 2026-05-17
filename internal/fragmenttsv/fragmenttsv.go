package fragmenttsv

import (
	"bufio"
	"fmt"
	"os"

	"github.com/ericksamera/radigest/internal/digest"
)

// Writer emits per-fragment TSV rows for downstream modeling. A Writer created
// with an empty path is a no-op, which lets callers keep TSV output disabled
// without nil checks.
type Writer struct {
	f         *os.File
	bw        *bufio.Writer
	closeFile bool
	disabled  bool
}

// New opens path and writes the TSV header. Use an empty path to disable TSV
// output. Use "-" to write to stdout.
func New(path string) (*Writer, error) {
	if path == "" {
		return &Writer{disabled: true}, nil
	}

	var f *os.File
	closeFile := false
	if path == "-" {
		f = os.Stdout
	} else {
		var err error
		f, err = os.Create(path)
		if err != nil {
			return nil, err
		}
		closeFile = true
	}
	w := &Writer{f: f, bw: bufio.NewWriter(f), closeFile: closeFile}
	if _, err := w.bw.WriteString("chrom\tstart0\tend0\tlength\thard_kept\tsize_weight\n"); err != nil {
		if closeFile {
			_ = f.Close()
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
	if w.closeFile {
		if closeErr := w.f.Close(); err == nil {
			err = closeErr
		}
	}
	return err
}
