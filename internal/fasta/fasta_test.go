package fasta

import (
	"os"
	"testing"
)

func TestStream(t *testing.T) {
	// tiny two-record FASTA (mixed case & spaces)
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
	for r := range ch {
		recs = append(recs, r)
	}
	if len(recs) != 2 {
		t.Fatalf("want 2 records, got %d", len(recs))
	}
	if recs[0].ID != "chr1" || string(recs[0].Seq) != "ACGTNN" {
		t.Fatalf("bad first record: %+v", recs[0])
	}
	if recs[1].ID != "chr2" || string(recs[1].Seq) != "GGCC" {
		t.Fatalf("bad second record: %+v", recs[1])
	}
}
