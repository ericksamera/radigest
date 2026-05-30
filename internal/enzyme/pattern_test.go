package enzyme

import "testing"

func TestStripCaret(t *testing.T) {
	site, cut := StripCaret("G^AATTC")
	if site != "GAATTC" || cut != 1 {
		t.Fatalf("with caret: got (%q,%d)", site, cut)
	}
	site, cut = StripCaret("AAAA")
	if site != "AAAA" || cut != 2 {
		t.Fatalf("no caret: got (%q,%d)", site, cut)
	}
}

func TestMatchMask_DegenerateR(t *testing.T) {
	mask := CompileMask("ACGTR") // R = A|G
	if !MatchMask(mask, []byte("ACGTA")) {
		t.Fatal("R should match A")
	}
	if !MatchMask(mask, []byte("ACGTG")) {
		t.Fatal("R should match G")
	}
	if MatchMask(mask, []byte("ACGTC")) {
		t.Fatal("R should not match C")
	}
}

func TestMatchMask_SequenceNDoesNotMatch(t *testing.T) {
	mask := CompileMask("A")
	if MatchMask(mask, []byte("N")) {
		t.Fatal("sequence 'N' must not match any site base")
	}
}

func TestCompileMaskCheckedRejectsInvalidIUPAC(t *testing.T) {
	if _, err := CompileMaskChecked("AX"); err == nil {
		t.Fatalf("CompileMaskChecked returned nil error for invalid IUPAC symbol")
	}
}

func TestIsExactACGT(t *testing.T) {
	if !IsExactACGT("GAATTC") {
		t.Fatalf("GAATTC should be recognized as an exact A/C/G/T motif")
	}
	if !IsExactACGT("gaattc") {
		t.Fatalf("lowercase exact motifs should be recognized")
	}
	if IsExactACGT("GCWGC") {
		t.Fatalf("degenerate motifs must not use the exact-match scanner")
	}
	if IsExactACGT("") {
		t.Fatalf("empty motif must not be exact")
	}
}

func TestBestMaskAnchorPrefersMostSpecificPosition(t *testing.T) {
	mask := CompileMask("NNACN")
	if got := BestMaskAnchor(mask); got != 2 {
		t.Fatalf("best anchor = %d, want 2", got)
	}
}

func TestMatchMaskAtUsesAnchorAndReferenceNStillFails(t *testing.T) {
	mask := CompileMask("NNACN")
	if !MatchMaskAt(mask, BestMaskAnchor(mask), []byte("GGACA")) {
		t.Fatalf("expected degenerate motif to match compatible reference")
	}
	if MatchMaskAt(mask, BestMaskAnchor(mask), []byte("GGNCG")) {
		t.Fatalf("reference N must not match even under degenerate motif N")
	}
}
