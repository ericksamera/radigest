package main

import (
	"math"
	"testing"
)

func TestValidatePositiveThreads(t *testing.T) {
	for _, n := range []int{1, 2, 16} {
		if err := validatePositiveThreads(n); err != nil {
			t.Fatalf("validatePositiveThreads(%d) returned error: %v", n, err)
		}
	}
	for _, n := range []int{0, -1} {
		if err := validatePositiveThreads(n); err == nil {
			t.Fatalf("validatePositiveThreads(%d) returned nil error", n)
		}
	}
}

func TestValidateSimGC(t *testing.T) {
	for _, gc := range []float64{0, 0.5, 1} {
		if err := validateSimGC(gc); err != nil {
			t.Fatalf("validateSimGC(%g) returned error: %v", gc, err)
		}
	}
	for _, gc := range []float64{-0.1, 1.1, math.NaN(), math.Inf(1), math.Inf(-1)} {
		if err := validateSimGC(gc); err == nil {
			t.Fatalf("validateSimGC(%g) returned nil error", gc)
		}
	}
}

func TestParseEnzymes(t *testing.T) {
	ens, names, err := parseEnzymes(" EcoRI , MseI ")
	if err != nil {
		t.Fatalf("parseEnzymes returned error: %v", err)
	}
	if len(ens) != 2 || ens[0].Name != "EcoRI" || ens[1].Name != "MseI" {
		t.Fatalf("unexpected enzymes: %+v", ens)
	}
	if len(names) != 2 || names[0] != "EcoRI" || names[1] != "MseI" {
		t.Fatalf("unexpected names: %+v", names)
	}
}

func TestParseEnzymesRejectsInvalidInputs(t *testing.T) {
	for _, value := range []string{
		"EcoRI,MseI,NcoI",
		"EcoRI,",
		"EcoRI,,MseI",
		"EcoRI,EcoRI",
		"NotAnEnzyme",
	} {
		if _, _, err := parseEnzymes(value); err == nil {
			t.Fatalf("parseEnzymes(%q) returned nil error", value)
		}
	}
}

func TestValidateOutputSelection(t *testing.T) {
	if err := validateOutputSelection("", "", "", "", "-"); err != nil {
		t.Fatalf("json stdout should satisfy output selection: %v", err)
	}
	if err := validateOutputSelection("", "", "", "", ""); err == nil {
		t.Fatalf("expected all disabled outputs to be rejected")
	}
}

func TestValidateOutputPathsRejectsCollisions(t *testing.T) {
	if err := validateOutputPaths("ref.fa", "ref.fa", "fragments.bed", "fragments.tsv", "fragments.fa", "run.json", true); err == nil {
		t.Fatalf("expected -gff/-fasta path collision to be rejected")
	}
	if err := validateOutputPaths("ref.fa", "out.gff3", "out.gff3", "fragments.tsv", "fragments.fa", "run.json", true); err == nil {
		t.Fatalf("expected -gff/-bed path collision to be rejected")
	}
	if err := validateOutputPaths("ref.fa", "out.gff3", "fragments.bed", "fragments.bed", "fragments.fa", "run.json", true); err == nil {
		t.Fatalf("expected -bed/-fragments-tsv path collision to be rejected")
	}
	if err := validateOutputPaths("ref.fa", "out.gff3", "fragments.bed", "fragments.tsv", "out.gff3", "run.json", true); err == nil {
		t.Fatalf("expected -gff/-fragments-fasta path collision to be rejected")
	}
	if err := validateOutputPaths("ref.fa", "out.gff3", "fragments.bed", "fragments.tsv", "ref.fa", "run.json", true); err == nil {
		t.Fatalf("expected -fasta/-fragments-fasta path collision to be rejected")
	}
	if err := validateOutputPaths("ref.fa", "out.gff3", "fragments.bed", "fragments.tsv", "fragments.fa", "out.gff3", true); err == nil {
		t.Fatalf("expected -gff/-json path collision to be rejected")
	}
	if err := validateOutputPaths("ref.fa", "-", "-", "fragments.tsv", "fragments.fa", "run.json", true); err == nil {
		t.Fatalf("expected two stdout outputs to be rejected")
	}
	if err := validateOutputPaths("ref.fa", "-", "", "", "fragments.fa", "-", true); err == nil {
		t.Fatalf("expected -gff stdout and -json stdout collision to be rejected")
	}
}

func TestValidateOutputPathsAllowsDisabledOutputsAndStdout(t *testing.T) {
	if err := validateOutputPaths("-", "-", "", "", "", "", true); err != nil {
		t.Fatalf("stdout/disabled outputs should be allowed: %v", err)
	}
	if err := validateOutputPaths("ref.fa", "", "", "", "", "-", true); err != nil {
		t.Fatalf("json stdout should be allowed: %v", err)
	}
	if err := validateOutputPaths("", "/dev/null", "", "", "", "", false); err != nil {
		t.Fatalf("/dev/null should be allowed as a sink: %v", err)
	}
}
