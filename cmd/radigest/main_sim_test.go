package main

import (
	"bytes"
	"encoding/json"
	"errors"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/ericksamera/radigest/internal/enzyme"
)

func TestSimFlags_RunWritesGFFJSONAndFragmentsTSV(t *testing.T) {
	dir := t.TempDir()
	gffPath := filepath.Join(dir, "frag.gff3")
	jsPath := filepath.Join(dir, "run.json")
	tsvPath := filepath.Join(dir, "fragments.tsv")

	var stdout, stderr bytes.Buffer
	if err := run([]string{
		"-sim-len", "10000",
		"-sim-gc", "0.50",
		"-sim-seed", "42",
		"-enzymes", "MluCI", // ^AATT (frequent; but we don't assert counts)
		"-gff", gffPath,
		"-fragments-tsv", tsvPath,
		"-json", jsPath,
		"-threads", "1",
	}, strings.NewReader(""), &stdout, &stderr); err != nil {
		t.Fatalf("run returned error: %v\nstderr:\n%s", err, stderr.String())
	}

	// GFF exists and has header
	gff, err := os.ReadFile(gffPath)
	if err != nil {
		t.Fatalf("read gff: %v", err)
	}
	if !bytes.HasPrefix(gff, []byte("##gff-version 3\n")) {
		t.Fatalf("gff missing header")
	}

	// TSV exists and has header
	tsv, err := os.ReadFile(tsvPath)
	if err != nil {
		t.Fatalf("read tsv: %v", err)
	}
	if !strings.HasPrefix(string(tsv), "chrom\tstart0\tend0\tlength\thard_kept\tsize_weight\n") {
		t.Fatalf("tsv missing header")
	}

	// JSON exists and is parseable with required fields
	var doc struct {
		Enzymes []string `json:"enzymes"`
		Min     int      `json:"min_length"`
		Max     int      `json:"max_length"`
		Size    struct {
			Model string `json:"model"`
		} `json:"size_selection"`
	}
	raw, err := os.ReadFile(jsPath)
	if err != nil {
		t.Fatalf("read json: %v", err)
	}
	if err := json.Unmarshal(raw, &doc); err != nil {
		t.Fatalf("parse json: %v", err)
	}
	if len(doc.Enzymes) == 0 || doc.Enzymes[0] != "MluCI" {
		t.Fatalf("json enzymes wrong: %+v", doc.Enzymes)
	}
	if doc.Size.Model != "hard" {
		t.Fatalf("json size model wrong: %q", doc.Size.Model)
	}
}

func TestRunRejectsInvalidEnzymeMetadataAsUsageError(t *testing.T) {
	old, existed := enzyme.DB["BadMeta"]
	enzyme.DB["BadMeta"] = enzyme.Enzyme{Name: "BadMeta", Recognition: "A^X", CutIndex: 1}
	defer func() {
		if existed {
			enzyme.DB["BadMeta"] = old
		} else {
			delete(enzyme.DB, "BadMeta")
		}
	}()

	var stdout, stderr bytes.Buffer
	err := run([]string{
		"-sim-len", "1000",
		"-sim-gc", "0.50",
		"-sim-seed", "42",
		"-enzymes", "BadMeta",
		"-threads", "1",
	}, strings.NewReader(""), &stdout, &stderr)
	if err == nil {
		t.Fatalf("run returned nil error for invalid enzyme metadata")
	}
	var usage usageError
	if !errors.As(err, &usage) {
		t.Fatalf("expected usageError, got %T: %v", err, err)
	}
	if exitCode(err) != 2 {
		t.Fatalf("expected exit code 2, got %d", exitCode(err))
	}
	if !strings.Contains(err.Error(), "invalid IUPAC") {
		t.Fatalf("expected invalid IUPAC error, got %v", err)
	}
}
