package main

import (
	"bytes"
	"encoding/json"
	"os"
	"path/filepath"
	"testing"
)

func TestRunCachedScreenWritesRankablePairJSON(t *testing.T) {
	dir := t.TempDir()
	fastaPath := filepath.Join(dir, "toy.fa")
	if err := os.WriteFile(fastaPath, []byte(">ecori_msei_double\nAAAAGAATTCTTAAAGAATTCTTT\n"), 0o644); err != nil {
		t.Fatalf("write FASTA: %v", err)
	}
	outDir := filepath.Join(dir, "screen")

	var stdout, stderr bytes.Buffer
	err := run([]string{
		"--fasta", fastaPath,
		"--enzymes", "EcoRI,MseI",
		"--min", "1",
		"--max", "100",
		"--score-min", "1",
		"--score-max", "100",
		"--size-model", "hard",
		"--out-dir", outDir,
		"--jobs", "1",
	}, &stdout, &stderr)
	if err != nil {
		t.Fatalf("run() error = %v\nstderr:\n%s", err, stderr.String())
	}

	jsonPath := filepath.Join(outDir, "json", "EcoRI__MseI.json")
	data, err := os.ReadFile(jsonPath)
	if err != nil {
		t.Fatalf("read pair JSON: %v", err)
	}
	var doc struct {
		Enzymes        []string `json:"enzymes"`
		TotalFragments int      `json:"total_fragments"`
		TotalBases     int      `json:"total_bases"`
		SizeSelection  struct {
			RawFragmentsInWindow int     `json:"raw_fragments_in_window"`
			RawBasesInWindow     int64   `json:"raw_bases_in_window"`
			WeightedFragments    float64 `json:"weighted_fragments"`
			WeightedBases        float64 `json:"weighted_bases"`
		} `json:"size_selection"`
		Screening struct {
			Engine string `json:"engine"`
		} `json:"screening"`
	}
	if err := json.Unmarshal(data, &doc); err != nil {
		t.Fatalf("unmarshal pair JSON: %v", err)
	}
	if got := doc.Enzymes; len(got) != 2 || got[0] != "EcoRI" || got[1] != "MseI" {
		t.Fatalf("enzymes = %#v, want EcoRI,MseI", got)
	}
	if doc.TotalFragments != 2 || doc.TotalBases != 11 {
		t.Fatalf("totals = fragments %d bases %d, want 2 and 11", doc.TotalFragments, doc.TotalBases)
	}
	if doc.SizeSelection.RawFragmentsInWindow != 2 || doc.SizeSelection.RawBasesInWindow != 11 {
		t.Fatalf("size-selection raw fields = %d/%d, want 2/11", doc.SizeSelection.RawFragmentsInWindow, doc.SizeSelection.RawBasesInWindow)
	}
	if doc.SizeSelection.WeightedFragments != 2 || doc.SizeSelection.WeightedBases != 11 {
		t.Fatalf("size-selection weighted fields = %g/%g, want 2/11", doc.SizeSelection.WeightedFragments, doc.SizeSelection.WeightedBases)
	}
	if doc.Screening.Engine != "cached-cut-index" {
		t.Fatalf("screening engine = %q, want cached-cut-index", doc.Screening.Engine)
	}
}

func TestDryRunDoesNotReadMissingFASTA(t *testing.T) {
	var stdout, stderr bytes.Buffer
	err := run([]string{
		"--fasta", filepath.Join(t.TempDir(), "missing.fa"),
		"--enzymes", "EcoRI,MseI",
		"--out-dir", t.TempDir(),
		"--dry-run",
	}, &stdout, &stderr)
	if err != nil {
		t.Fatalf("dry-run returned error: %v", err)
	}
	if !bytes.Contains(stdout.Bytes(), []byte("EcoRI__MseI")) {
		t.Fatalf("dry-run stdout did not contain pair tag; stdout=%q", stdout.String())
	}
}

func TestWriteJSONAtomicWritesCompleteFile(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "pair.json")

	if err := writeJSONAtomic(path, map[string]any{"ok": true}); err != nil {
		t.Fatalf("writeJSONAtomic returned error: %v", err)
	}

	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("read output JSON: %v", err)
	}
	var doc map[string]bool
	if err := json.Unmarshal(data, &doc); err != nil {
		t.Fatalf("unmarshal output JSON: %v", err)
	}
	if !doc["ok"] {
		t.Fatalf("decoded JSON = %#v, want ok=true", doc)
	}

	matches, err := filepath.Glob(filepath.Join(dir, "pair.json.*.tmp"))
	if err != nil {
		t.Fatalf("glob temp files: %v", err)
	}
	if len(matches) != 0 {
		t.Fatalf("temporary files remain after successful write: %#v", matches)
	}
}

func TestWriteJSONAtomicDoesNotLeaveFinalFileOnEncodeError(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "pair.json")

	err := writeJSONAtomic(path, map[string]any{"bad": make(chan int)})
	if err == nil {
		t.Fatalf("writeJSONAtomic returned nil error for unencodable value")
	}
	if _, statErr := os.Stat(path); !os.IsNotExist(statErr) {
		t.Fatalf("final JSON file exists after failed write; stat error = %v", statErr)
	}

	matches, globErr := filepath.Glob(filepath.Join(dir, "pair.json.*.tmp"))
	if globErr != nil {
		t.Fatalf("glob temp files: %v", globErr)
	}
	if len(matches) != 0 {
		t.Fatalf("temporary files remain after failed write: %#v", matches)
	}
}
