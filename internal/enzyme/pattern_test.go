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
	if !MatchMask(mask, []byte("ACGTA")) { t.Fatal("R should match A") }
	if !MatchMask(mask, []byte("ACGTG")) { t.Fatal("R should match G") }
	if  MatchMask(mask, []byte("ACGTC")) { t.Fatal("R should not match C") }
}

func TestMatchMask_SequenceNDoesNotMatch(t *testing.T) {
	mask := CompileMask("A")
	if MatchMask(mask, []byte("N")) {
		t.Fatal("sequence 'N' must not match any site base")
	}
}
