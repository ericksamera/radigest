package main

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strings"

	"github.com/ericksamera/radigest/internal/enzyme"
)

func validatePositiveThreads(n int) error {
	if n < 1 {
		return fmt.Errorf("-threads must be >= 1 (got %d)", n)
	}
	return nil
}

func validateSimGC(gc float64) error {
	if math.IsNaN(gc) || math.IsInf(gc, 0) || gc < 0 || gc > 1 {
		return fmt.Errorf("-sim-gc must be in [0,1] (got %g)", gc)
	}
	return nil
}

func parseEnzymes(value string) ([]enzyme.Enzyme, []string, error) {
	parts := strings.Split(value, ",")
	names := make([]string, 0, len(parts))
	for _, part := range parts {
		name := strings.TrimSpace(part)
		if name == "" {
			return nil, nil, fmt.Errorf("invalid -enzymes %q: empty enzyme name", value)
		}
		names = append(names, name)
	}
	if len(names) > 2 {
		return nil, nil, fmt.Errorf("invalid -enzymes %q: specify one or two enzymes", value)
	}

	ens := make([]enzyme.Enzyme, 0, len(names))
	canonicalNames := make([]string, 0, len(names))
	for _, name := range names {
		e, ok := enzyme.DB[name]
		if !ok {
			return nil, nil, fmt.Errorf("unknown enzyme %q", name)
		}
		ens = append(ens, e)
		canonicalNames = append(canonicalNames, e.Name)
	}
	if len(ens) == 2 && ens[0].Name == ens[1].Name {
		return nil, nil, fmt.Errorf("first two enzymes must differ (got %s,%s)", ens[0].Name, ens[1].Name)
	}
	return ens, canonicalNames, nil
}

func validateOutputSelection(gffPath, bedPath, fragmentsTSVPath, fragmentsFASTAPath, jsonPath string) error {
	for _, path := range []string{gffPath, bedPath, fragmentsTSVPath, fragmentsFASTAPath, jsonPath} {
		if activeOutputPath(path) {
			return nil
		}
	}
	return fmt.Errorf("no outputs enabled; omit output flags to write JSON summary to stdout, or set -json, -gff, -bed, -fragments-tsv, or -fragments-fasta")
}

func validateOutputPaths(fastaPath, gffPath, bedPath, fragmentsTSVPath, fragmentsFASTAPath, jsonPath string, hasFastaInput bool) error {
	outputs := []namedPath{
		{name: "-gff", path: gffPath, stdoutAllowed: true},
		{name: "-bed", path: bedPath, stdoutAllowed: true},
		{name: "-fragments-tsv", path: fragmentsTSVPath, stdoutAllowed: true},
		{name: "-fragments-fasta", path: fragmentsFASTAPath, stdoutAllowed: true},
		{name: "-json", path: jsonPath, stdoutAllowed: true},
	}

	if hasFastaInput && activeFilePath(fastaPath) {
		inputKey, err := comparablePath(fastaPath)
		if err != nil {
			return fmt.Errorf("-fasta path: %w", err)
		}
		for _, out := range outputs {
			if !activeOutputPath(out.path) {
				continue
			}
			outKey, err := outputComparablePath(out)
			if err != nil {
				return fmt.Errorf("%s path: %w", out.name, err)
			}
			if inputKey == outKey {
				return fmt.Errorf("refusing to use input FASTA %q as output for %s", fastaPath, out.name)
			}
		}
	}

	seen := make(map[string]string)
	for _, out := range outputs {
		if !activeOutputPath(out.path) {
			continue
		}
		key, err := outputComparablePath(out)
		if err != nil {
			return fmt.Errorf("%s path: %w", out.name, err)
		}
		if prior, ok := seen[key]; ok {
			return fmt.Errorf("refusing to write both %s and %s to %q", prior, out.name, out.path)
		}
		seen[key] = out.name
	}
	return nil
}

type namedPath struct {
	name          string
	path          string
	stdoutAllowed bool
}

func activeFilePath(path string) bool {
	return strings.TrimSpace(path) != "" && path != "-"
}

func activeOutputPath(path string) bool {
	return strings.TrimSpace(path) != ""
}

func outputComparablePath(path namedPath) (string, error) {
	if path.path == "-" && path.stdoutAllowed {
		return "<stdout>", nil
	}
	return comparablePath(path.path)
}

func comparablePath(path string) (string, error) {
	abs, err := filepath.Abs(path)
	if err != nil {
		return "", err
	}
	abs = filepath.Clean(abs)
	if resolved, err := filepath.EvalSymlinks(abs); err == nil {
		abs = resolved
	} else if !os.IsNotExist(err) {
		return "", err
	}
	return abs, nil
}
