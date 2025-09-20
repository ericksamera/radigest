package fasta

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"io"
	"os"
)

const bufSize = 4 << 20 // 4 MiB

// Record is one FASTA entry (whole chromosome or contig).
type Record struct {
	ID  string
	Seq []byte // upper-case, no newlines; reused – copy if you need to keep it
}

// Stream reads `path` (file path or "-" for STDIN) and sends each record.
// It closes the channel when done or on first error (returned).
func Stream(path string, out chan<- Record) error {
	var src io.Reader

	if path == "-" {
		br := bufio.NewReader(os.Stdin)
		if magic, _ := br.Peek(2); len(magic) == 2 && magic[0] == 0x1f && magic[1] == 0x8b {
			gz, err := gzip.NewReader(br)
			if err != nil {
				return err
			}
			defer gz.Close()
			src = gz
		} else {
			src = br
		}
	} else {
		f, err := os.Open(path)
		if err != nil {
			return err
		}
		defer f.Close()
		br := bufio.NewReader(f)
		if magic, _ := br.Peek(2); len(magic) == 2 && magic[0] == 0x1f && magic[1] == 0x8b {
			gz, err := gzip.NewReader(br)
			if err != nil {
				return err
			}
			defer gz.Close()
			src = gz
		} else {
			src = br
		}
	}

	r := bufio.NewReaderSize(src, bufSize)
	var (
		id   []byte
		seq  []byte
	)
	flush := func() {
		if id != nil {
			out <- Record{ID: string(id), Seq: bytes.ToUpper(seq)}
			seq = seq[:0]
		}
	}
	for {
		line, err := r.ReadBytes('\n')
		if err != nil && err != io.EOF {
			return err
		}
		if len(line) > 0 && line[len(line)-1] == '\n' {
			line = line[:len(line)-1] // trim newline
		}
		if len(line) > 0 && line[0] == '>' { // header
			flush()
			id = bytes.Fields(line[1:])[0] // up to first space
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
