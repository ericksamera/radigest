package fasta

import (
	"bytes"
	"compress/gzip"
	"os"
	"testing"
)

func TestStream_StdinPlain(t *testing.T) {
	r, w, _ := os.Pipe()
	old := os.Stdin
	os.Stdin = r
	defer func() { os.Stdin = old }()

	go func() {
		_, _ = w.Write([]byte(">chr1\nACGT\n>chr2\nNN\n"))
		_ = w.Close()
	}()

	ch := make(chan Record)
	go func() { if err := Stream("-", ch); err != nil { t.Error(err) } }()

	var recs []Record
	for r := range ch { recs = append(recs, r) }
	if len(recs) != 2 || recs[0].ID != "chr1" || string(recs[1].Seq) != "NN" {
		t.Fatalf("bad records: %+v", recs)
	}
}

func TestStream_StdinGzip(t *testing.T) {
	var buf bytes.Buffer
	zw := gzip.NewWriter(&buf)
	_, _ = zw.Write([]byte(">chr1\nac\n>chr2\nGg\n"))
	_ = zw.Close()

	r, w, _ := os.Pipe()
	old := os.Stdin
	os.Stdin = r
	defer func() { os.Stdin = old }()

	go func() { _, _ = w.Write(buf.Bytes()); _ = w.Close() }()

	ch := make(chan Record)
	go func() { if err := Stream("-", ch); err != nil { t.Error(err) } }()
	var recs []Record
	for r := range ch { recs = append(recs, r) }
	if len(recs) != 2 || string(recs[0].Seq) != "AC" || string(recs[1].Seq) != "GG" {
		t.Fatalf("bad gunzip via stdin: %+v", recs)
	}
}

func TestStream_WindowsCRLF_Trimmed(t *testing.T) {
	data := []byte(">c\r\nA\r\nC\r\n")
	tmp, _ := os.CreateTemp("", "crlf*.fa")
	defer os.Remove(tmp.Name())
	_, _ = tmp.Write(data); _ = tmp.Close()

	ch := make(chan Record)
	go func() { if err := Stream(tmp.Name(), ch); err != nil { t.Error(err) } }()
	rec := <-ch
	if bytes.Contains(rec.Seq, []byte{'\r'}) {
		t.Fatal("CR should be trimmed from sequences")
	}
	for range ch {}
}
