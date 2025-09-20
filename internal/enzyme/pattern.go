package enzyme

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
func CompileMask(site string) []uint8 {
	b := []byte(site)
	m := make([]uint8, len(b))
	for i, c := range b {
		if c >= 'a' && c <= 'z' { // upper-case on the fly
			c -= 'a' - 'A'
		}
		m[i] = codeMap[c]
	}
	return m
}

// baseMaskWin maps a reference base to its mask for matching.
// NOTE: We *block* 'N' in the reference (mask=0) so 'N' never matches any site.
func baseMaskWin(b byte) uint8 {
	if b >= 'a' && b <= 'z' { b -= 'a' - 'A' }
	if b == 'N' {
		return 0
	}
	if m, ok := codeMap[b]; ok {
		return m
	}
	return 0 // anything unknown in the sequence fails to match
}

// MatchMask returns true iff window matches the compiled mask.
func MatchMask(mask []uint8, window []byte) bool {
	n := len(mask)
	// fast reject on last position
	if baseMaskWin(window[n-1])&mask[n-1] == 0 {
		return false
	}
	for i := 0; i < n-1; i++ {
		if baseMaskWin(window[i])&mask[i] == 0 {
			return false
		}
	}
	return true
}
