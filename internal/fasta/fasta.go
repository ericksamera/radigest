package fasta

import (
	"bufio"
	"bytes"
	// "fmt"
	"io"
	"os"
)

const bufSize = 4 << 20 // 4 MiB

// Record is one FASTA entry (whole chromosome or contig).
type Record struct {
	ID   string
	Seq  []byte // upper-case, no newlines; reused – copy if you need to keep it
}

// Stream reads `path` and sends each record down the chan.
// It closes the channel when done or on first error (returned).
func Stream(path string, out chan<- Record) error {
	f, err := os.Open(path)
	if err != nil {
		return err
	}
	defer f.Close()

	r := bufio.NewReaderSize(f, bufSize)
	var (
		id   []byte
		seq  []byte
		line []byte
	)
	flush := func() {
		if id != nil {
			out <- Record{ID: string(id), Seq: bytes.ToUpper(seq)}
			seq = seq[:0]
		}
	}
	for {
		line, err = r.ReadBytes('\n')
		if err != nil && err != io.EOF {
			return err
		}
		if len(line) > 0 && line[len(line)-1] == '\n' {
			line = line[:len(line)-1] // trim newline
		}
		if len(line) > 0 && line[0] == '>' { // header
			flush()
			id = bytes.Fields(line[1:])[0] // grab up-to-first-space
		} else {
			seq = append(seq, line...)
		}
		if err == io.EOF {
			flush()
			close(out)
			return nil
		}
	}
}

// Example usage (in cmd/main.go):
//
// ch := make(chan fasta.Record, 2)
// go func() {
//     if err := fasta.Stream(fastaPath, ch); err != nil {
//         log.Fatal(err)
//     }
// }()
// for rec := range ch {
//     // hand rec.Seq to digest workers
// }
