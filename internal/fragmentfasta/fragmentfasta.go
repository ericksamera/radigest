package fragmentfasta

import (
	"bufio"
	"fmt"
	"io"
	"os"

	"github.com/ericksamera/radigest/internal/digest"
	"github.com/ericksamera/radigest/internal/gff"
)

const wrapWidth = 80

// Writer emits FASTA records for fragments. A Writer created with an empty path
// is a no-op, which lets callers keep FASTA output disabled without nil checks.
type Writer struct {
	bw       *bufio.Writer
	close    func() error
	disabled bool
}

// New opens path for fragment FASTA output. Use an empty path to disable FASTA
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
	return &Writer{bw: bufio.NewWriter(sink), close: close}, nil
}

// Write emits one fragment FASTA record. Coordinates are 0-based half-open in
// the header, and ordinal should match the corresponding saved fragment ordinal
// for the chromosome.
func (w *Writer) Write(chr string, ordinal int, fr digest.Fragment, seq []byte) error {
	if w == nil || w.disabled {
		return nil
	}
	if fr.Start < 0 || fr.End < fr.Start || fr.End > len(seq) {
		return fmt.Errorf(
			"fragment FASTA: invalid fragment for %s_%d: start=%d end=%d sequence_length=%d",
			chr,
			ordinal,
			fr.Start,
			fr.End,
			len(seq),
		)
	}

	length := fr.End - fr.Start
	escapedChr := gff.EscapeAttributeValue(chr)
	if escapedChr == "" {
		escapedChr = "."
	}
	if _, err := fmt.Fprintf(
		w.bw,
		">%s chrom=%s start0=%d end0=%d length=%d\n",
		fragmentID(chr, ordinal),
		escapedChr,
		fr.Start,
		fr.End,
		length,
	); err != nil {
		return err
	}

	fragmentSeq := seq[fr.Start:fr.End]
	for len(fragmentSeq) > 0 {
		n := wrapWidth
		if len(fragmentSeq) < n {
			n = len(fragmentSeq)
		}
		if _, err := w.bw.Write(fragmentSeq[:n]); err != nil {
			return err
		}
		if err := w.bw.WriteByte('\n'); err != nil {
			return err
		}
		fragmentSeq = fragmentSeq[n:]
	}
	return nil
}

// Close flushes pending FASTA output and closes owned files. Stdout is flushed
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

func fragmentID(chr string, ordinal int) string {
	if chr == "" {
		return fmt.Sprintf("frag%d", ordinal)
	}
	return fmt.Sprintf("%s_%d", gff.EscapeAttributeValue(chr), ordinal)
}
