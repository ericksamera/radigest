package collector

import (
	"bufio"
	"os"
	"strings"
	"testing"

	"github.com/KPU-AGC/radigest/internal/digest"
)

func TestWriterWriteStream(t *testing.T) {
	tmp, err := os.CreateTemp("", "frag*.gff")
	if err != nil {
		t.Fatal(err)
	}
	defer os.Remove(tmp.Name())
	_ = tmp.Close()

	w, err := NewWriter(tmp.Name())
	if err != nil {
		t.Fatal(err)
	}
	frags := make(chan digest.Fragment, 2)
	frags <- digest.Fragment{Start: 0, End: 5}
	frags <- digest.Fragment{Start: 10, End: 12}
	close(frags)

	cs, err := w.WriteStream("chr1", frags)
	if err != nil {
		t.Fatal(err)
	}
	stats, err := w.Close()
	if err != nil {
		t.Fatal(err)
	}
	if cs.Fragments != 2 || cs.Bases != 7 || stats.TotalFragments != 2 || stats.TotalBases != 7 {
		t.Fatalf("stats wrong: chr=%+v total=%+v", cs, stats)
	}

	f, err := os.Open(tmp.Name())
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	sc := bufio.NewScanner(f)
	lines := []string{}
	for sc.Scan() {
		lines = append(lines, sc.Text())
	}
	if err := sc.Err(); err != nil {
		t.Fatal(err)
	}
	if len(lines) != 3 {
		t.Fatalf("expected header + 2 lines, got %d", len(lines))
	}
	if !strings.Contains(lines[1], "ID=chr1_1;Length=5") || !strings.Contains(lines[2], "ID=chr1_2;Length=2") {
		t.Fatalf("unexpected GFF lines:\n%s\n%s", lines[1], lines[2])
	}
}
