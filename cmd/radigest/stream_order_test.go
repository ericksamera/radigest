package main

import (
	"bufio"
	"os"
	"strings"
	"testing"

	"github.com/KPU-AGC/radigest/internal/collector"
	"github.com/KPU-AGC/radigest/internal/digest"
)

func result(idx int, chr string, frags ...digest.Fragment) digestResult {
	fragCh := make(chan digest.Fragment, len(frags))
	for _, fr := range frags {
		fragCh <- fr
	}
	close(fragCh)
	errCh := make(chan error, 1)
	errCh <- nil
	close(errCh)
	return digestResult{idx: idx, chr: chr, frags: fragCh, errors: errCh}
}

func TestWriteResultStreamsWritesChromosomesInIndexOrder(t *testing.T) {
	tmp, err := os.CreateTemp("", "frag*.gff")
	if err != nil {
		t.Fatal(err)
	}
	defer func() { _ = os.Remove(tmp.Name()) }()
	_ = tmp.Close()

	w, err := collector.NewWriter(tmp.Name())
	if err != nil {
		t.Fatal(err)
	}
	results := make(chan digestResult, 2)
	results <- result(1, "chr2", digest.Fragment{Start: 10, End: 12})
	results <- result(0, "chr1", digest.Fragment{Start: 0, End: 5})
	close(results)

	if err := writeResultStreams(w, results, false); err != nil {
		t.Fatal(err)
	}
	stats, err := w.Close()
	if err != nil {
		t.Fatal(err)
	}
	if stats.TotalFragments != 2 || stats.TotalBases != 7 {
		t.Fatalf("bad stats: %+v", stats)
	}

	f, err := os.Open(tmp.Name())
	if err != nil {
		t.Fatal(err)
	}
	defer func() { _ = f.Close() }()
	sc := bufio.NewScanner(f)
	lines := []string{}
	for sc.Scan() {
		lines = append(lines, sc.Text())
	}
	if err := sc.Err(); err != nil {
		t.Fatal(err)
	}
	if len(lines) < 3 {
		t.Fatalf("expected header + 2 lines, got %d", len(lines))
	}
	if !strings.HasPrefix(lines[1], "chr1\t") || !strings.HasPrefix(lines[2], "chr2\t") {
		t.Fatalf("order wrong:\n%s\n%s", lines[1], lines[2])
	}
}
