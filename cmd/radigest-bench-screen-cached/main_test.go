package main

import (
	"bytes"
	"encoding/csv"
	"os"
	"path/filepath"
	"testing"
)

func TestRunWritesPhaseTimingTSV(t *testing.T) {
	dir := t.TempDir()
	fastaPath := filepath.Join(dir, "toy.fa")
	if err := os.WriteFile(fastaPath, []byte(">ecori_msei_double\nAAAAGAATTCTTAAAGAATTCTTT\n"), 0o644); err != nil {
		t.Fatalf("write FASTA: %v", err)
	}

	var stdout, stderr bytes.Buffer
	err := run([]string{
		"--fasta", fastaPath,
		"--enzymes", "EcoRI,MseI,PstI",
		"--min", "1",
		"--max", "100",
		"--score-min", "1",
		"--score-max", "100",
		"--size-model", "hard",
		"--jobs", "2",
		"--runs", "1",
		"--output-mode", "marshal",
	}, &stdout, &stderr)
	if err != nil {
		t.Fatalf("run() error = %v\nstderr:\n%s", err, stderr.String())
	}

	reader := csv.NewReader(bytes.NewReader(stdout.Bytes()))
	reader.Comma = '\t'
	rows, err := reader.ReadAll()
	if err != nil {
		t.Fatalf("parse TSV: %v\nstdout:\n%s", err, stdout.String())
	}
	if len(rows) != 2 {
		t.Fatalf("got %d TSV rows, want header+1; stdout:\n%s", len(rows), stdout.String())
	}
	header := indexHeader(rows[0])
	if rows[1][header["candidate_enzymes"]] != "3" {
		t.Fatalf("candidate_enzymes = %s, want 3", rows[1][header["candidate_enzymes"]])
	}
	if rows[1][header["candidate_pairs"]] != "3" {
		t.Fatalf("candidate_pairs = %s, want 3", rows[1][header["candidate_pairs"]])
	}
	if rows[1][header["summaries"]] != "3" {
		t.Fatalf("summaries = %s, want 3", rows[1][header["summaries"]])
	}
	if rows[1][header["output_mode"]] != "marshal" {
		t.Fatalf("output_mode = %s, want marshal", rows[1][header["output_mode"]])
	}
	if rows[1][header["digest_checksum"]] == "" {
		t.Fatalf("digest_checksum is empty")
	}
}

func TestRunWriteModeWritesJSONFiles(t *testing.T) {
	dir := t.TempDir()
	fastaPath := filepath.Join(dir, "toy.fa")
	outDir := filepath.Join(dir, "bench")
	if err := os.WriteFile(fastaPath, []byte(">ecori_msei_double\nAAAAGAATTCTTAAAGAATTCTTT\n"), 0o644); err != nil {
		t.Fatalf("write FASTA: %v", err)
	}

	var stdout, stderr bytes.Buffer
	err := run([]string{
		"--fasta", fastaPath,
		"--enzymes", "EcoRI,MseI",
		"--min", "1",
		"--max", "100",
		"--score-min", "1",
		"--score-max", "100",
		"--size-model", "hard",
		"--jobs", "1",
		"--runs", "1",
		"--output-mode", "write",
		"--out-dir", outDir,
	}, &stdout, &stderr)
	if err != nil {
		t.Fatalf("run() error = %v\nstderr:\n%s", err, stderr.String())
	}

	path := filepath.Join(outDir, "run_001", "json", "EcoRI__MseI.json")
	if _, err := os.Stat(path); err != nil {
		t.Fatalf("expected output JSON at %s: %v", path, err)
	}
}

func indexHeader(header []string) map[string]int {
	out := make(map[string]int, len(header))
	for i, name := range header {
		out[name] = i
	}
	return out
}
