// internal/enzyme/iupac.go
package enzyme

import "fmt"

// 4‑bit mask per base
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
    'N': 0, // ← was 1<<0|1<<1|1<<2|1<<3
}

// CompilePattern converts an IUPAC recognition string to a slice of 4‑bit masks.
func CompilePattern(seq string) ([]uint8, error) {
    out := make([]uint8, len(seq))
    for i := 0; i < len(seq); i++ {
        c := seq[i]
        if c >= 'a' && c <= 'z' {         // ensure upper-case
            c -= 'a' - 'A'
        }
        m, ok := codeMap[c]
        if !ok {
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
		base := window[i]
        if base >= 'a' && base <= 'z' {
            base -= 'a' - 'A'
        }
        if (m & codeMap[base]) == 0 {
            return false
        }
    }
    return true
}
