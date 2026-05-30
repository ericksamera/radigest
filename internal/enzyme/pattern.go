package enzyme

import "math/bits"

// StripCaret removes “^” from the recognition site and returns (cleanSite, cutOffset).
func StripCaret(recog string) (string, int) {
	for i := 0; i < len(recog); i++ {
		if recog[i] == '^' {
			return recog[:i] + recog[i+1:], i
		}
	}
	// no caret: default cut mid-site
	return recog, len(recog) / 2
}

// CompileMask converts an IUPAC string to per-position bit-masks.
//
// CompileMask preserves the historical unchecked behavior: unknown symbols are
// converted to zero masks. Use CompileMaskChecked for user- or database-facing
// validation.
func CompileMask(site string) []uint8 {
	b := []byte(site)
	m := make([]uint8, len(b))
	for i, c := range b {
		m[i] = motifMaskTable[c]
	}
	return m
}

// CompileMaskChecked converts an IUPAC string to per-position bit-masks and
// rejects unknown symbols instead of silently producing zero masks.
func CompileMaskChecked(site string) ([]uint8, error) {
	return CompilePattern(site)
}

// baseMaskWin maps a reference base to its mask for matching.
// NOTE: We *block* 'N' in the reference (mask=0) so 'N' never matches any site.
func baseMaskWin(b byte) uint8 {
	return refMaskTable[b]
}

// BestMaskAnchor returns the most selective position in a compiled motif mask.
// Positions with fewer allowed bases are tested first by MatchMaskAt.
func BestMaskAnchor(mask []uint8) int {
	if len(mask) == 0 {
		return 0
	}
	best := 0
	bestPop := 9
	for i, m := range mask {
		pop := bits.OnesCount8(m)
		if pop < bestPop {
			best = i
			bestPop = pop
			if pop == 1 {
				break
			}
		}
	}
	return best
}

// MatchMask returns true iff window matches the compiled mask.
func MatchMask(mask []uint8, window []byte) bool {
	return MatchMaskAt(mask, BestMaskAnchor(mask), window)
}

// MatchMaskAt returns true iff window matches the compiled mask, testing the
// provided anchor position first for a fast reject.
func MatchMaskAt(mask []uint8, anchor int, window []byte) bool {
	n := len(mask)
	if n == 0 || len(window) < n {
		return false
	}
	if anchor < 0 || anchor >= n {
		anchor = n - 1
	}
	if baseMaskWin(window[anchor])&mask[anchor] == 0 {
		return false
	}
	for i := 0; i < n; i++ {
		if i == anchor {
			continue
		}
		if baseMaskWin(window[i])&mask[i] == 0 {
			return false
		}
	}
	return true
}

// IsExactACGT reports whether site contains only unambiguous A/C/G/T bases.
func IsExactACGT(site string) bool {
	if site == "" {
		return false
	}
	for i := 0; i < len(site); i++ {
		switch site[i] {
		case 'A', 'C', 'G', 'T', 'a', 'c', 'g', 't':
			continue
		default:
			return false
		}
	}
	return true
}
