package main

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestMainWritesFragmentFASTAForHardKeptFragments(t *testing.T) {
	dir := t.TempDir()
	refPath := filepath.Join(dir, "ref.fa")
	gffPath := filepath.Join(dir, "frag.gff3")
	fastaPath := filepath.Join(dir, "fragments.fa")
	if err := os.WriteFile(refPath, []byte(">chr1\nAAAAGAATTCTTAAAGAATTC\n"), 0o644); err != nil {
		t.Fatal(err)
	}

	var stdout, stderr bytes.Buffer
	if err := run([]string{
		"-fasta", refPath,
		"-enzymes", "EcoRI,MseI",
		"-min", "1",
		"-max", "1000",
		"-gff", gffPath,
		"-fragments-tsv", "",
		"-fragments-fasta", fastaPath,
		"-threads", "1",
	}, strings.NewReader(""), &stdout, &stderr); err != nil {
		t.Fatalf("run returned error: %v\nstderr:\n%s", err, stderr.String())
	}

	raw, err := os.ReadFile(fastaPath)
	if err != nil {
		t.Fatalf("read fragment FASTA: %v", err)
	}
	got := string(raw)
	want := strings.Join([]string{
		">chr1_1 chrom=chr1 start0=5 end0=11 length=6",
		"AATTCT",
		">chr1_2 chrom=chr1 start0=11 end0=16 length=5",
		"TAAAG",
		"",
	}, "\n")
	if got != want {
		t.Fatalf("fragment FASTA mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}
