package fasta

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"fmt"
	"io"
	"os"
)

const bufSize = 4 << 20 // 4 MiB

type Record struct {
	ID  string
	Seq []byte // upper-case, no newlines; reused – copy if you need to keep it
}

// Stream reads path (file path or "-" for STDIN) and sends each record.
// It always closes out before returning, including when it returns an error.
func Stream(path string, out chan<- Record) (err error) {
	if out == nil {
		return fmt.Errorf("fasta stream: output channel is nil")
	}
	defer close(out)

	src, cleanup, err := openSource(path)
	if err != nil {
		return err
	}
	defer func() {
		if closeErr := cleanup(); err == nil && closeErr != nil {
			err = closeErr
		}
	}()

	r := bufio.NewReaderSize(src, bufSize)
	var id, seq []byte
	seenHeader := false
	lineNo := 0

	flush := func() {
		if id != nil {
			out <- Record{ID: string(id), Seq: bytes.ToUpper(seq)}
			seq = seq[:0]
		}
	}

	for {
		line, readErr := r.ReadBytes('\n')
		if readErr != nil && readErr != io.EOF {
			return fmt.Errorf("fasta stream: read %q: %w", path, readErr)
		}
		if len(line) > 0 {
			lineNo++
			line = bytes.TrimRight(line, "\r\n")

			switch {
			case len(line) == 0:
				// Ignore blank lines. They are common enough in hand-edited FASTA
				// files and do not change sequence content.
			case line[0] == '>':
				flush()
				fields := bytes.Fields(line[1:])
				if len(fields) == 0 {
					return fmt.Errorf("fasta stream: empty header at line %d", lineNo)
				}
				id = append([]byte(nil), fields[0]...)
				seenHeader = true
			case id == nil:
				return fmt.Errorf("fasta stream: sequence before first header at line %d", lineNo)
			default:
				seq = append(seq, line...)
			}
		}
		if readErr == io.EOF {
			if !seenHeader {
				return fmt.Errorf("fasta stream: no records in %q", path)
			}
			flush()
			return nil
		}
	}
}

func openSource(path string) (io.Reader, func() error, error) {
	if path == "-" {
		br := bufio.NewReader(os.Stdin)
		if isGzip(br) {
			gz, err := gzip.NewReader(br)
			if err != nil {
				return nil, nil, err
			}
			return gz, gz.Close, nil
		}
		return br, func() error { return nil }, nil
	}

	f, err := os.Open(path)
	if err != nil {
		return nil, nil, err
	}
	br := bufio.NewReader(f)
	if isGzip(br) {
		gz, err := gzip.NewReader(br)
		if err != nil {
			_ = f.Close()
			return nil, nil, err
		}
		return gz, func() error {
			gzErr := gz.Close()
			fileErr := f.Close()
			if gzErr != nil {
				return gzErr
			}
			return fileErr
		}, nil
	}
	return br, f.Close, nil
}

func isGzip(br *bufio.Reader) bool {
	magic, _ := br.Peek(2)
	return len(magic) == 2 && magic[0] == 0x1f && magic[1] == 0x8b
}
