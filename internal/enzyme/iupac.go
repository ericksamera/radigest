package enzyme

import "fmt"

// 4-bit mask per base.
var codeMap = map[byte]uint8{
	'A': 1 << 0,
	'C': 1 << 1,
	'G': 1 << 2,
	'T': 1 << 3,
	'R': (1 << 0) | (1 << 2),
	'Y': (1 << 1) | (1 << 3),
	'S': (1 << 1) | (1 << 2),
	'W': (1 << 0) | (1 << 3),
	'K': (1 << 2) | (1 << 3),
	'M': (1 << 0) | (1 << 1),
	'B': (1 << 1) | (1 << 2) | (1 << 3),
	'D': (1 << 0) | (1 << 2) | (1 << 3),
	'H': (1 << 0) | (1 << 1) | (1 << 3),
	'V': (1 << 0) | (1 << 1) | (1 << 2),
	'N': (1 << 0) | (1 << 1) | (1 << 2) | (1 << 3), // match any base
}

// motifMaskTable maps IUPAC motif symbols to their allowed A/C/G/T masks.
// Motif N means any unambiguous A/C/G/T reference base.
var motifMaskTable [256]uint8

// refMaskTable maps reference bases to masks used while matching a reference
// window. Reference N and all unrecognized symbols intentionally map to zero,
// so no recognition site is inferred across assembly gaps or unknown bases.
var refMaskTable [256]uint8

func init() {
	for b, m := range codeMap {
		setMaskBothCases(&motifMaskTable, b, m)
	}

	setMaskBothCases(&refMaskTable, 'A', 1<<0)
	setMaskBothCases(&refMaskTable, 'C', 1<<1)
	setMaskBothCases(&refMaskTable, 'G', 1<<2)
	setMaskBothCases(&refMaskTable, 'T', 1<<3)

	// Reference N remains zero by design.
	refMaskTable['N'] = 0
	refMaskTable['n'] = 0
}

func setMaskBothCases(table *[256]uint8, b byte, mask uint8) {
	table[b] = mask
	if b >= 'A' && b <= 'Z' {
		table[b+'a'-'A'] = mask
	}
}

// CompilePattern converts an IUPAC recognition string to a slice of 4-bit masks.
func CompilePattern(seq string) ([]uint8, error) {
	out := make([]uint8, len(seq))
	for i := 0; i < len(seq); i++ {
		c := seq[i]
		m := motifMaskTable[c]
		if m == 0 {
			return nil, fmt.Errorf("invalid IUPAC base %q", c)
		}
		out[i] = m
	}
	return out, nil
}

// Match tests whether window matches pattern (case-sensitive, same length).
func Match(pattern []uint8, window []byte) bool {
	if len(pattern) != len(window) {
		return false
	}
	for i, m := range pattern {
		if m&refMaskTable[window[i]] == 0 {
			return false
		}
	}
	return true
}
