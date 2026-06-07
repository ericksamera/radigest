package main

import (
	"bytes"
	"encoding/json"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestNoOutputFlagsWritesJSONToStdoutOnly(t *testing.T) {
	dir := t.TempDir()
	refPath := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(refPath, []byte(">chr1\nAAAAGAATTCTTAAAGAATTC\n"), 0o644); err != nil {
		t.Fatal(err)
	}

	oldDir, err := os.Getwd()
	if err != nil {
		t.Fatal(err)
	}
	if err := os.Chdir(dir); err != nil {
		t.Fatal(err)
	}
	defer func() { _ = os.Chdir(oldDir) }()

	stdout, _ := runCaptured(t, []string{
		"-fasta", refPath,
		"-enzymes", "EcoRI,MseI",
		"-min", "1",
		"-max", "1000",
		"-threads", "1",
	}, "")

	var doc struct {
		Enzymes        []string `json:"enzymes"`
		TotalFragments int      `json:"total_fragments"`
		TotalBases     int      `json:"total_bases"`
	}
	if err := json.Unmarshal([]byte(stdout), &doc); err != nil {
		t.Fatalf("parse stdout JSON: %v\nstdout:\n%s", err, stdout)
	}
	if len(doc.Enzymes) != 2 || doc.Enzymes[0] != "EcoRI" || doc.Enzymes[1] != "MseI" {
		t.Fatalf("json enzymes wrong: %+v", doc.Enzymes)
	}
	if doc.TotalFragments != 2 || doc.TotalBases != 11 {
		t.Fatalf("json stats wrong: %+v", doc)
	}
	for _, unexpected := range []string{"fragments.gff3", "fragments.bed", "fragments.tsv"} {
		if _, err := os.Stat(filepath.Join(dir, unexpected)); !os.IsNotExist(err) {
			t.Fatalf("unexpected default output file %s exists or stat failed: %v", unexpected, err)
		}
	}
}

func TestExplicitJSONStdoutDoesNotCreateDefaultArtifacts(t *testing.T) {
	dir := t.TempDir()
	refPath := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(refPath, []byte(">chr1\nAAAAGAATTCTTAAAGAATTC\n"), 0o644); err != nil {
		t.Fatal(err)
	}

	oldDir, err := os.Getwd()
	if err != nil {
		t.Fatal(err)
	}
	if err := os.Chdir(dir); err != nil {
		t.Fatal(err)
	}
	defer func() { _ = os.Chdir(oldDir) }()

	stdout, _ := runCaptured(t, []string{
		"-fasta", refPath,
		"-enzymes", "EcoRI,MseI",
		"-json", "-",
		"-threads", "1",
	}, "")

	if !json.Valid([]byte(stdout)) {
		t.Fatalf("stdout is not valid JSON: %q", stdout)
	}
	for _, unexpected := range []string{"fragments.gff3", "fragments.bed", "fragments.tsv"} {
		if _, err := os.Stat(filepath.Join(dir, unexpected)); !os.IsNotExist(err) {
			t.Fatalf("unexpected default output file %s exists or stat failed: %v", unexpected, err)
		}
	}
}

func TestSummaryIncludesSchemaProvenanceAndResolvedSimSeed(t *testing.T) {
	args := []string{
		"-sim-len", "10000",
		"-sim-gc", "0.50",
		"-sim-seed", "0",
		"-enzymes", "MluCI",
		"-threads", "1",
	}
	stdout, _ := runCaptured(t, args, "")

	var doc struct {
		SchemaVersion   int      `json:"schema_version"`
		RadigestVersion string   `json:"radigest_version"`
		Command         []string `json:"command"`
		Input           struct {
			Source           string   `json:"source"`
			SimLength        int      `json:"sim_length"`
			SimGC            *float64 `json:"sim_gc"`
			SimSeedRequested *int64   `json:"sim_seed_requested"`
			SimSeedResolved  *int64   `json:"sim_seed_resolved"`
		} `json:"input"`
		Parameters struct {
			MinLength   int    `json:"min_length"`
			MaxLength   int    `json:"max_length"`
			ScoreMin    int    `json:"score_min"`
			ScoreMax    int    `json:"score_max"`
			SizeModel   string `json:"size_model"`
			Threads     int    `json:"threads"`
			AllowSame   bool   `json:"allow_same"`
			StrictCuts  bool   `json:"strict_cuts"`
			IncludeEnds bool   `json:"include_ends"`
		} `json:"parameters"`
		Outputs struct {
			JSON string `json:"json"`
		} `json:"outputs"`
		Warnings []string `json:"warnings"`
	}
	if err := json.Unmarshal([]byte(stdout), &doc); err != nil {
		t.Fatalf("parse stdout JSON: %v\nstdout:\n%s", err, stdout)
	}
	if doc.SchemaVersion != summarySchemaVersion {
		t.Fatalf("schema_version: got %d want %d", doc.SchemaVersion, summarySchemaVersion)
	}
	if doc.RadigestVersion != version {
		t.Fatalf("radigest_version: got %q want %q", doc.RadigestVersion, version)
	}
	wantCommand := append([]string{"radigest"}, args...)
	if strings.Join(doc.Command, "\x00") != strings.Join(wantCommand, "\x00") {
		t.Fatalf("command mismatch\nwant: %+v\ngot:  %+v", wantCommand, doc.Command)
	}
	if doc.Input.Source != "simulation" || doc.Input.SimLength != 10000 {
		t.Fatalf("input summary wrong: %+v", doc.Input)
	}
	if doc.Input.SimGC == nil || *doc.Input.SimGC != 0.50 {
		t.Fatalf("sim_gc not recorded correctly: %+v", doc.Input.SimGC)
	}
	if doc.Input.SimSeedRequested == nil || *doc.Input.SimSeedRequested != 0 {
		t.Fatalf("requested sim seed not recorded correctly: %+v", doc.Input.SimSeedRequested)
	}
	if doc.Input.SimSeedResolved == nil || *doc.Input.SimSeedResolved == 0 {
		t.Fatalf("resolved sim seed not recorded correctly: %+v", doc.Input.SimSeedResolved)
	}
	if doc.Parameters.MinLength != 1 || doc.Parameters.MaxLength != 1<<30 || doc.Parameters.ScoreMin != 1 || doc.Parameters.ScoreMax != 1<<30 {
		t.Fatalf("length parameters wrong: %+v", doc.Parameters)
	}
	if doc.Parameters.SizeModel != "hard" || doc.Parameters.Threads != 1 || doc.Parameters.AllowSame || doc.Parameters.StrictCuts || doc.Parameters.IncludeEnds {
		t.Fatalf("parameters wrong: %+v", doc.Parameters)
	}
	if doc.Outputs.JSON != "-" {
		t.Fatalf("default JSON output not recorded as stdout: %+v", doc.Outputs)
	}
	if !containsWarning(doc.Warnings, "time-based seed") {
		t.Fatalf("warnings missing time-based seed notice: %+v", doc.Warnings)
	}
}

func TestRunWritesBEDForHardKeptFragmentsAndRecordsOutput(t *testing.T) {
	dir := t.TempDir()
	refPath := filepath.Join(dir, "ref.fa")
	bedPath := filepath.Join(dir, "fragments.bed")
	jsonPath := filepath.Join(dir, "run.json")
	if err := os.WriteFile(refPath, []byte(">chr1\nAAAAGAATTCTTAAAGAATTC\n"), 0o644); err != nil {
		t.Fatal(err)
	}

	var stdout, stderr bytes.Buffer
	if err := run([]string{
		"-fasta", refPath,
		"-enzymes", "EcoRI,MseI",
		"-min", "1",
		"-max", "1000",
		"-bed", bedPath,
		"-json", jsonPath,
		"-threads", "1",
	}, strings.NewReader(""), &stdout, &stderr); err != nil {
		t.Fatalf("run returned error: %v\nstderr:\n%s", err, stderr.String())
	}

	raw, err := os.ReadFile(bedPath)
	if err != nil {
		t.Fatalf("read BED: %v", err)
	}
	want := strings.Join([]string{
		"chr1\t5\t11\tchr1_1\t0\t+",
		"chr1\t11\t16\tchr1_2\t0\t+",
		"",
	}, "\n")
	if string(raw) != want {
		t.Fatalf("BED mismatch\nwant:\n%s\ngot:\n%s", want, string(raw))
	}

	data, err := os.ReadFile(jsonPath)
	if err != nil {
		t.Fatalf("read JSON: %v", err)
	}
	var doc struct {
		BED     string `json:"bed"`
		Outputs struct {
			BED  string `json:"bed"`
			JSON string `json:"json"`
		} `json:"outputs"`
	}
	if err := json.Unmarshal(data, &doc); err != nil {
		t.Fatalf("parse JSON: %v\njson:\n%s", err, string(data))
	}
	if doc.BED != bedPath || doc.Outputs.BED != bedPath || doc.Outputs.JSON != jsonPath {
		t.Fatalf("BED/JSON output paths not recorded correctly: %+v", doc)
	}
}

func TestRunReadsInjectedStdinAndWritesInjectedStdout(t *testing.T) {
	stdout, _ := runCaptured(t, []string{
		"-fasta", "-",
		"-enzymes", "EcoRI,MseI",
		"-min", "1",
		"-max", "1000",
		"-gff", "-",
		"-threads", "1",
	}, ">chr1\nAAAAGAATTCTTAAAGAATTC\n")

	if !strings.HasPrefix(stdout, "##gff-version 3\n") {
		t.Fatalf("stdout missing GFF header: %q", stdout)
	}
	if !strings.Contains(stdout, "chr1\tradigest\tfragment\t6\t11") {
		t.Fatalf("stdout missing expected fragment: %q", stdout)
	}
}

func TestRunMissingRequiredFlagsReturnsUsageError(t *testing.T) {
	var stdout, stderr bytes.Buffer
	err := run([]string{"-fasta", "ref.fa"}, strings.NewReader(""), &stdout, &stderr)
	if err == nil {
		t.Fatalf("expected error")
	}
	if exitCode(err) != 2 {
		t.Fatalf("expected usage exit code 2, got %d for %v", exitCode(err), err)
	}
	if !strings.Contains(stderr.String(), "Required flags:") {
		t.Fatalf("stderr missing usage text: %q", stderr.String())
	}
}

func runCaptured(t *testing.T, args []string, stdin string) (string, string) {
	t.Helper()

	var stdout, stderr bytes.Buffer
	if err := run(args, strings.NewReader(stdin), &stdout, &stderr); err != nil {
		t.Fatalf("run returned error: %v\nstderr:\n%s", err, stderr.String())
	}
	return stdout.String(), stderr.String()
}

func containsWarning(warnings []string, needle string) bool {
	for _, warning := range warnings {
		if strings.Contains(warning, needle) {
			return true
		}
	}
	return false
}
