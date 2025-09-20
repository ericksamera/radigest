package main

import (
	"bytes"
	"encoding/json"
	"flag"
	"os"
	"path/filepath"
	"testing"
)

func TestSimFlags_RunWritesGFFAndJSON(t *testing.T) {
	dir := t.TempDir()
	gffPath := filepath.Join(dir, "frag.gff3")
	jsPath  := filepath.Join(dir, "run.json")

	// Fresh FlagSet so main() can define flags
	flag.CommandLine = flag.NewFlagSet("radigest", flag.ExitOnError)
	os.Args = []string{
		"radigest",
		"-sim-len", "10000",
		"-sim-gc", "0.50",
		"-sim-seed", "42",
		"-enzymes", "MluCI", // ^AATT (frequent; but we don't assert counts)
		"-gff", gffPath,
		"-json", jsPath,
		"-threads", "1",
	}
	main()

	// GFF exists and has header
	gff, err := os.ReadFile(gffPath)
	if err != nil { t.Fatalf("read gff: %v", err) }
	if !bytes.HasPrefix(gff, []byte("##gff-version 3\n")) {
		t.Fatalf("gff missing header")
	}

	// JSON exists and is parseable with required fields
	var doc struct{
		Enzymes []string `json:"enzymes"`
		Min     int      `json:"min_length"`
		Max     int      `json:"max_length"`
	}
	raw, err := os.ReadFile(jsPath)
	if err != nil { t.Fatalf("read json: %v", err) }
	if err := json.Unmarshal(raw, &doc); err != nil {
		t.Fatalf("parse json: %v", err)
	}
	if len(doc.Enzymes) == 0 || doc.Enzymes[0] != "MluCI" {
		t.Fatalf("json enzymes wrong: %+v", doc.Enzymes)
	}
}
