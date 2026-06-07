package main

import (
	"bytes"
	"encoding/csv"
	"encoding/json"
	"errors"
	"os"
	"path/filepath"
	"strings"
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
	summaryTSVPath := filepath.Join(outDir, "design.summary.tsv")
	jsonPath := filepath.Join(outDir, "design.json")
	reportPath := filepath.Join(outDir, "design.report.txt")
	rows := readTSV(t, tsvPath)
	if len(rows) != 4 {
		t.Fatalf("got %d TSV rows, want header + 3 candidates", len(rows))
	}
	header := indexHeader(rows[0])
	for _, name := range []string{"fit_score", "fit_loss", "predicted_weighted_genome_pct", "target_mean_locus_depth", "predicted_mean_locus_depth"} {
		if _, ok := header[name]; !ok {
			t.Fatalf("TSV header missing %s: %#v", name, rows[0])
		}
	}
	for _, oldName := range []string{"design_score", "design_loss", "generated_weighted_genome_pct", "desired_depth", "expected_mean_depth"} {
		if _, ok := header[oldName]; ok {
			t.Fatalf("TSV header still contains old field %s: %#v", oldName, rows[0])
		}
	}
	if rows[1][header["enzyme_a"]] != "EcoRI" || rows[1][header["enzyme_b"]] != "MseI" {
		t.Fatalf("top pair = %s,%s, want EcoRI,MseI", rows[1][header["enzyme_a"]], rows[1][header["enzyme_b"]])
	}
	if rows[1][header["feasible"]] != "true" {
		t.Fatalf("top pair feasible = %s, want true", rows[1][header["feasible"]])
	}

	summaryRows := readTSV(t, summaryTSVPath)
	if len(summaryRows) != 4 {
		t.Fatalf("got %d summary TSV rows, want header + 3 candidates", len(summaryRows))
	}
	summaryHeader := indexHeader(summaryRows[0])
	for _, name := range []string{"rank", "enzyme_pair", "feasible", "target_pct", "predicted_pct", "target_depth", "predicted_depth", "max_samples", "fit_score"} {
		if _, ok := summaryHeader[name]; !ok {
			t.Fatalf("summary TSV header missing %s: %#v", name, summaryRows[0])
		}
	}
	if summaryRows[1][summaryHeader["enzyme_pair"]] != "EcoRI,MseI" {
		t.Fatalf("summary top enzyme_pair = %s, want EcoRI,MseI", summaryRows[1][summaryHeader["enzyme_pair"]])
	}
	if summaryRows[1][summaryHeader["feasible"]] != "true" {
		t.Fatalf("summary top feasible = %s, want true", summaryRows[1][summaryHeader["feasible"]])
	}

	reportRows := readKeyValueReport(t, reportPath)
	for key, want := range map[string]string{
		"format":        "radigest-design-report-v1",
		"best_pair":     "EcoRI,MseI",
		"best_feasible": "true",
		"output.json":   jsonPath,
		"output.report": reportPath,
	} {
		if got := reportRows[key]; got != want {
			t.Fatalf("report[%s] = %q, want %q; report=%+v", key, got, want, reportRows)
		}
	}
	if reportRows["best_predicted_weighted_genome_pct"] == "" || reportRows["best_predicted_mean_locus_depth"] == "" {
		t.Fatalf("report missing predicted design metrics: %+v", reportRows)
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
		Outputs struct {
			TSV        string `json:"tsv"`
			SummaryTSV string `json:"summary_tsv"`
			JSON       string `json:"json"`
			Report     string `json:"report"`
		} `json:"outputs"`
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
	if report.Outputs.TSV != tsvPath || report.Outputs.SummaryTSV != summaryTSVPath || report.Outputs.JSON != jsonPath || report.Outputs.Report != reportPath {
		t.Fatalf("output paths not recorded correctly: %+v", report.Outputs)
	}
	for _, needle := range []string{"size_model\thard", "hard_size_window_bp\t1-100", "score_range_bp\t1-100", "size_mean_bp\tNA", "size_sd_bp\tNA"} {
		if !strings.Contains(stderr.String(), needle) {
			t.Fatalf("stderr missing %q; stderr:\n%s", needle, stderr.String())
		}
	}
	if !strings.Contains(stderr.String(), "design_summary_tsv\t"+summaryTSVPath) {
		t.Fatalf("stderr missing design summary path %q; stderr:\n%s", summaryTSVPath, stderr.String())
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
	if _, err := os.Stat(filepath.Join(outDir, "design.summary.tsv")); err != nil {
		t.Fatalf("expected design.summary.tsv: %v", err)
	}
	if _, err := os.Stat(filepath.Join(outDir, "design.json")); err != nil {
		t.Fatalf("expected design.json: %v", err)
	}
	if _, err := os.Stat(filepath.Join(outDir, "design.report.txt")); err != nil {
		t.Fatalf("expected design.report.txt: %v", err)
	}
}

func TestRunHelpShowsGroupedDesignHelp(t *testing.T) {
	var stdout, stderr bytes.Buffer
	err := run([]string{"--help"}, &stdout, &stderr)
	if err != nil {
		t.Fatalf("help returned error: %v", err)
	}
	help := stderr.String()
	for _, needle := range []string{
		"radigest-design\n",
		"Author:  Gennerick J. Samera (erick.samera@kpu.ca)",
		"Version: " + version,
		"License: MIT",
		"Description: rank enzyme pairs from genome-fraction, depth, and sequencing-budget targets",
		"Required design inputs:",
		"Size selection and recovery model:",
		"Sequencing budget:",
		"Ranking and scoring:",
		"Size-selection models:",
		"normal",
		"triangular",
		"soft-window",
		"Outputs written:",
		"design.report.txt",
	} {
		if !strings.Contains(help, needle) {
			t.Fatalf("help missing %q; help:\n%s", needle, help)
		}
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

func readKeyValueReport(t *testing.T, path string) map[string]string {
	t.Helper()
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("read report: %v", err)
	}
	out := make(map[string]string)
	for _, line := range strings.Split(strings.TrimSpace(string(data)), "\n") {
		parts := strings.SplitN(line, "\t", 2)
		if len(parts) != 2 {
			t.Fatalf("invalid report line %q in report:\n%s", line, string(data))
		}
		out[parts[0]] = parts[1]
	}
	return out
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
