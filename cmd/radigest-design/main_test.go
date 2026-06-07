package main

import (
	"bytes"
	"encoding/csv"
	"encoding/json"
	"errors"
	"os"
	"path/filepath"
	"testing"
)

func TestRunWritesDesignOutputs(t *testing.T) {
	dir := t.TempDir()
	fastaPath := filepath.Join(dir, "toy.fa")
	if err := os.WriteFile(fastaPath, []byte(">ecori_msei_double\nAAAAGAATTCTTAAAGAATTCTTT\n"), 0o644); err != nil {
		t.Fatalf("write FASTA: %v", err)
	}
	outDir := filepath.Join(dir, "design")

	var stdout, stderr bytes.Buffer
	err := run([]string{
		"--fasta", fastaPath,
		"--enzymes", "EcoRI,MseI,PstI",
		"--min", "1",
		"--max", "100",
		"--score-min", "1",
		"--score-max", "100",
		"--size-model", "hard",
		"--target-genome-pct", "45.833333",
		"--coverage-tolerance-pct", "1",
		"--desired-depth", "10",
		"--samples", "1",
		"--read-layout", "pe",
		"--read-length", "150",
		"--lane-read-pairs", "1000",
		"--out-dir", outDir,
		"--jobs", "1",
		"--build-workers", "2",
	}, &stdout, &stderr)
	if err != nil {
		t.Fatalf("run() error = %v\nstderr:\n%s", err, stderr.String())
	}

	tsvPath := filepath.Join(outDir, "design.tsv")
	jsonPath := filepath.Join(outDir, "design.json")
	rows := readTSV(t, tsvPath)
	if len(rows) != 4 {
		t.Fatalf("got %d TSV rows, want header + 3 candidates", len(rows))
	}
	header := indexHeader(rows[0])
	if rows[1][header["enzyme_a"]] != "EcoRI" || rows[1][header["enzyme_b"]] != "MseI" {
		t.Fatalf("top pair = %s,%s, want EcoRI,MseI", rows[1][header["enzyme_a"]], rows[1][header["enzyme_b"]])
	}
	if rows[1][header["feasible"]] != "true" {
		t.Fatalf("top pair feasible = %s, want true", rows[1][header["feasible"]])
	}

	data, err := os.ReadFile(jsonPath)
	if err != nil {
		t.Fatalf("read JSON: %v", err)
	}
	var report struct {
		SchemaVersion int      `json:"schema_version"`
		Command       []string `json:"command"`
		Summary       struct {
			CandidatePairs int      `json:"candidate_pairs"`
			FeasiblePairs  int      `json:"feasible_pairs"`
			BestPair       []string `json:"best_pair"`
		} `json:"summary"`
		Results []struct {
			EnzymeA  string `json:"enzyme_a"`
			EnzymeB  string `json:"enzyme_b"`
			Feasible bool   `json:"feasible"`
		} `json:"results"`
	}
	if err := json.Unmarshal(data, &report); err != nil {
		t.Fatalf("unmarshal JSON: %v\n%s", err, string(data))
	}
	if report.SchemaVersion != 1 || report.Summary.CandidatePairs != 3 || len(report.Results) != 3 {
		t.Fatalf("unexpected report summary: %+v", report)
	}
	if len(report.Command) == 0 || report.Command[0] != "radigest-design" {
		t.Fatalf("command[0] = %#v, want radigest-design", report.Command)
	}
	if got := report.Summary.BestPair; len(got) != 2 || got[0] != "EcoRI" || got[1] != "MseI" {
		t.Fatalf("best_pair = %+v, want EcoRI,MseI", got)
	}
}

func TestRunAcceptsConciseDesignAliases(t *testing.T) {
	dir := t.TempDir()
	fastaPath := filepath.Join(dir, "toy.fa")
	if err := os.WriteFile(fastaPath, []byte(">ecori_msei_double\nAAAAGAATTCTTAAAGAATTCTTT\n"), 0o644); err != nil {
		t.Fatalf("write FASTA: %v", err)
	}
	outDir := filepath.Join(dir, "design")

	var stdout, stderr bytes.Buffer
	err := run([]string{
		"--ref", fastaPath,
		"--enzymes", "EcoRI,MseI",
		"--min", "1",
		"--max", "100",
		"--score-min", "1",
		"--score-max", "100",
		"--size-model", "hard",
		"--pct", "45.833333",
		"--coverage-tolerance-pct", "1",
		"--depth", "10",
		"--samples", "1",
		"--read-layout", "pe",
		"--read-length", "150",
		"--flowcell-read-pairs", "1000",
		"--out-dir", outDir,
		"--jobs", "1",
		"--build-workers", "2",
	}, &stdout, &stderr)
	if err != nil {
		t.Fatalf("run() error = %v\nstderr:\n%s", err, stderr.String())
	}

	if _, err := os.Stat(filepath.Join(outDir, "design.tsv")); err != nil {
		t.Fatalf("expected design.tsv: %v", err)
	}
	if _, err := os.Stat(filepath.Join(outDir, "design.json")); err != nil {
		t.Fatalf("expected design.json: %v", err)
	}
}

func TestRunRejectsMissingBudget(t *testing.T) {
	var stdout, stderr bytes.Buffer
	err := run([]string{
		"--fasta", "ref.fa",
		"--enzymes", "EcoRI,MseI",
		"--target-genome-pct", "2.5",
		"--desired-depth", "10",
		"--samples", "1",
		"--read-length", "150",
	}, &stdout, &stderr)
	if err == nil {
		t.Fatalf("expected error")
	}
	if exitCode(err) != 2 {
		t.Fatalf("exitCode = %d, want 2", exitCode(err))
	}
}

func TestParsePositiveCount(t *testing.T) {
	got, err := parsePositiveCount("300M")
	if err != nil {
		t.Fatalf("parsePositiveCount returned error: %v", err)
	}
	if got != 300000000 {
		t.Fatalf("got %g, want 300000000", got)
	}
	gotInt, err := parsePositiveCountInt("2_643_888_753")
	if err != nil {
		t.Fatalf("parsePositiveCountInt returned error: %v", err)
	}
	if gotInt != 2643888753 {
		t.Fatalf("got %d, want 2643888753", gotInt)
	}
}

func readTSV(t *testing.T, path string) [][]string {
	t.Helper()
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("read TSV: %v", err)
	}
	reader := csv.NewReader(bytes.NewReader(data))
	reader.Comma = '\t'
	rows, err := reader.ReadAll()
	if err != nil {
		t.Fatalf("parse TSV: %v\n%s", err, string(data))
	}
	return rows
}

func indexHeader(header []string) map[string]int {
	out := make(map[string]int, len(header))
	for i, name := range header {
		out[name] = i
	}
	return out
}

func exitCode(err error) int {
	var usage usageError
	if errors.As(err, &usage) {
		return 2
	}
	return 1
}
