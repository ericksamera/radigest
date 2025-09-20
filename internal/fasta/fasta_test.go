package fasta

import (
	"bytes"
	"compress/gzip"
	"os"
	"testing"
)

func TestStream(t *testing.T) {
	data := ">chr1\nacgT\nNN\n>chr2 some desc\nGgCc\n"

	tmp, err := os.CreateTemp("", "fasta*.fa")
	if err != nil { t.Fatal(err) }
	defer os.Remove(tmp.Name())
	if _, err = tmp.WriteString(data); err != nil { t.Fatal(err) }
	tmp.Close()

	ch := make(chan Record)
	go func() {
		if err := Stream(tmp.Name(), ch); err != nil {
			t.Error(err)
		}
	}()

	var recs []Record
	for r := range ch { recs = append(recs, r) }
	if len(recs) != 2 { t.Fatalf("want 2 records, got %d", len(recs)) }
	if recs[0].ID != "chr1" || string(recs[0].Seq) != "ACGTNN" { t.Fatalf("bad chr1: %+v", recs[0]) }
	if recs[1].ID != "chr2" || string(recs[1].Seq) != "GGCC"   { t.Fatalf("bad chr2: %+v", recs[1]) }
}

func TestStreamGzip(t *testing.T) {
	data := ">chr1\nacgT\nNN\n>chr2 some desc\nGgCc\n"
	tmp, err := os.CreateTemp("", "fasta*.fa.gz")
	if err != nil { t.Fatal(err) }
	defer os.Remove(tmp.Name())

	var buf bytes.Buffer
	zw := gzip.NewWriter(&buf)
	if _, err := zw.Write([]byte(data)); err != nil { t.Fatal(err) }
	if err := zw.Close(); err != nil { t.Fatal(err) }
	if _, err := tmp.Write(buf.Bytes()); err != nil { t.Fatal(err) }
	tmp.Close()

	ch := make(chan Record)
	go func() {
		if err := Stream(tmp.Name(), ch); err != nil {
			t.Error(err)
		}
	}()
	var recs []Record
	for r := range ch { recs = append(recs, r) }
	if len(recs) != 2 { t.Fatalf("want 2 records, got %d", len(recs)) }
	if recs[0].ID != "chr1" || string(recs[0].Seq) != "ACGTNN" { t.Fatalf("bad chr1: %+v", recs[0]) }
	if recs[1].ID != "chr2" || string(recs[1].Seq) != "GGCC"   { t.Fatalf("bad chr2: %+v", recs[1]) }
}
