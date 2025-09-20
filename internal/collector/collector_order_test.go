package collector

import (
	"bufio"
	"os"
	"strings"
	"testing"

	"radigest/internal/digest"
)

func TestCollector_OutOfOrderIdxWritesInOrder(t *testing.T) {
	tmp, _ := os.CreateTemp("", "frag*.gff")
	defer os.Remove(tmp.Name())

	in, done, err := New(tmp.Name())
	if err != nil { t.Fatal(err) }

	// Send idx=1 first, then idx=0; collector must serialize 0 then 1.
	in <- Msg{Idx: 1, Chr: "chr2", Frags: []digest.Fragment{{Start: 10, End: 12}}}
	in <- Msg{Idx: 0, Chr: "chr1", Frags: []digest.Fragment{{Start: 0, End: 5}}}
	close(in)
	<-done

	f, _ := os.Open(tmp.Name())
	defer f.Close()
	sc := bufio.NewScanner(f)

	lines := []string{}
	for sc.Scan() { lines = append(lines, sc.Text()) }
	if len(lines) < 3 {
		t.Fatalf("expected header + 2 lines, got %d", len(lines))
	}
	if !strings.HasPrefix(lines[1], "chr1\t") || !strings.HasPrefix(lines[2], "chr2\t") {
		t.Fatalf("order wrong:\n%s\n%s", lines[1], lines[2])
	}
}

func TestCollector_EmptyThenNonEmpty(t *testing.T) {
	tmp, _ := os.CreateTemp("", "frag*.gff")
	defer os.Remove(tmp.Name())

	in, done, _ := New(tmp.Name())
	in <- Msg{Idx: 0, Chr: "empty", Frags: nil}
	in <- Msg{Idx: 1, Chr: "chrX", Frags: []digest.Fragment{{Start: 1, End: 4}}}
	close(in)
	stats := <-done

	if stats.TotalFragments != 1 || stats.PerChr["chrX"].Fragments != 1 {
		t.Fatalf("stats wrong: %+v", stats)
	}
}
