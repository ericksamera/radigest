package fragmenttsv

import (
	"bufio"
	"fmt"
	"os"

	"github.com/KPU-AGC/radigest/internal/digest"
)

type Writer struct {
	f         *os.File
	bw        *bufio.Writer
	closeFile bool
}

func New(path string) (*Writer, error) {
	if path == "" {
		return nil, nil
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

func (w *Writer) Write(chr string, fr digest.Fragment, hardKept bool, weight float64) error {
	if w == nil {
		return nil
	}
	length := fr.End - fr.Start
	_, err := fmt.Fprintf(w.bw, "%s\t%d\t%d\t%d\t%t\t%.8g\n", chr, fr.Start, fr.End, length, hardKept, weight)
	return err
}

func (w *Writer) Close() error {
	if w == nil {
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
