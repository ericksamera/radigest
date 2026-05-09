package main

import (
	"fmt"
	"math"
	"strings"

	"github.com/KPU-AGC/radigest/internal/enzyme"
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
