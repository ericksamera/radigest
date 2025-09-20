package collector

import (
	"os"
	"testing"

	"radigest/internal/digest"
)

func TestCollector(t *testing.T) {
	tmp, _ := os.CreateTemp("", "frag*.gff")
	defer os.Remove(tmp.Name())

	in, done, err := New(tmp.Name())
	if err != nil { t.Fatal(err) }

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
