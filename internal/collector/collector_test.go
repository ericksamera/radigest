package collector

import (
	"os"
	"strings"
	"testing"

	"github.com/ericksamera/radigest/internal/digest"
)

func TestCollector(t *testing.T) {
	tmp, err := os.CreateTemp("", "frag*.gff")
	if err != nil {
		t.Fatal(err)
	}
	defer func() { _ = os.Remove(tmp.Name()) }()
	if err := tmp.Close(); err != nil {
		t.Fatal(err)
	}

	in, done, err := New(tmp.Name())
	if err != nil {
		t.Fatal(err)
	}

	// send two chromosomes in deterministic order (Idx 0,1)
	in <- Msg{Idx: 0, Chr: "chr1", Frags: []digest.Fragment{{Start: 0, End: 5}}}
	in <- Msg{Idx: 1, Chr: "chr2", Frags: []digest.Fragment{
		{Start: 10, End: 15}, {Start: 20, End: 26},
	}}
	close(in)

	stats := <-done
	if stats.TotalFragments != 3 || stats.TotalBases != 16 {
		t.Fatalf("bad stats: %+v", stats)
	}
	if stats.PerChr["chr2"].Fragments != 2 {
		t.Fatalf("per-chr stats wrong: %+v", stats.PerChr)
	}
}

func TestWriterEscapesGFF3Fields(t *testing.T) {
	tmp, err := os.CreateTemp("", "frag*.gff")
	if err != nil {
		t.Fatal(err)
	}
	defer func() { _ = os.Remove(tmp.Name()) }()
	if err := tmp.Close(); err != nil {
		t.Fatal(err)
	}

	w, err := NewWriter(tmp.Name())
	if err != nil {
		t.Fatal(err)
	}
	chr := "chr 1;ID=x,50%\tbad"
	if err := w.WriteFragment(chr, 1, digest.Fragment{Start: 0, End: 5}); err != nil {
		t.Fatal(err)
	}
	if _, err := w.Close(); err != nil {
		t.Fatal(err)
	}

	raw, err := os.ReadFile(tmp.Name())
	if err != nil {
		t.Fatal(err)
	}
	got := string(raw)
	if !strings.Contains(got, "chr%201%3BID%3Dx%2C50%25%09bad\tradigest") {
		t.Fatalf("seqid was not escaped: %q", got)
	}
	if !strings.Contains(got, "ID=chr%201%3BID%3Dx%2C50%25%09bad_1;Length=5") {
		t.Fatalf("attributes were not escaped: %q", got)
	}
}
