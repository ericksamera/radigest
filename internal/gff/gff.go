package gff

import (
	"bufio"
	"fmt"
	"io"
	"os"

	"radigest/internal/digest"
)

// WriteFile writes one chromosome worth of fragments to `path`.
// Creates or truncates the file; caller passes chromosome ID.
func WriteFile(path, chr string, frags []digest.Fragment) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	return Write(f, chr, frags)
}

// Write streams a minimal, valid GFF3 (version header + one line / fragment).
// Coordinates are converted to *1-based closed* as GFF expects.
func Write(w io.Writer, chr string, frags []digest.Fragment) error {
	bw := bufio.NewWriter(w)
	if _, err := bw.WriteString("##gff-version 3\n"); err != nil {
		return err
	}
	for i, f := range frags {
		start := f.Start + 1 // 1-based
		end   := f.End       // half-open → closed
		if _, err := fmt.Fprintf(
			bw,
			"%s\tradigest\tfragment\t%d\t%d\t.\t+\t.\tID=frag%d\n",
			chr, start, end, i+1,
		); err != nil {
			return err
		}
	}
	return bw.Flush()
}
