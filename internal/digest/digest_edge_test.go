package digest

import (
	"testing"

	"radigest/internal/enzyme"
)

func TestSingleDigest_NoOrOneCutKeepsNone(t *testing.T) {
	eA := enzyme.DB["EcoRI"] // G^AATTC
	if got := Digest([]byte("TTTGGGTTT"), []enzyme.Enzyme{eA}, 0, 1<<30); len(got) != 0 {
		t.Fatalf("no site: got %d", len(got))
	}
	if got := Digest([]byte("TTTGAATTCTTT"), []enzyme.Enzyme{eA}, 0, 1<<30); len(got) != 0 {
		t.Fatalf("one site: want 0 frags, got %d", len(got))
	}
}

func TestSingleDigest_InclusiveMinMax(t *testing.T) {
	eA := enzyme.DB["MluCI"] // ^AATT
	seq := []byte("AATTCCCCAATT") // cuts at 0 and 8 → len 8
	if got := Digest(seq, []enzyme.Enzyme{eA}, 8, 8); len(got) != 1 {
		t.Fatalf("min=max=8 should keep exactly 1 frag, got %d", len(got))
	}
	if got := Digest(seq, []enzyme.Enzyme{eA}, 9, 9); len(got) != 0 {
		t.Fatalf("length 9 should keep none, got %d", len(got))
	}
}

func TestDoubleDigest_AAorBBSuppressed(t *testing.T) {
	eA := enzyme.DB["EcoRI"]
	eB := enzyme.DB["NcoI"] // absent
	seq := []byte("AAAAGAATTCAAAAGAATTCAAA")
	if got := Digest(seq, []enzyme.Enzyme{eA, eB}, 1, 1<<30); len(got) != 0 {
		t.Fatalf("AA adjacency must not be kept by default, got %d", len(got))
	}
}

func TestDoubleDigest_ZeroLengthWhenCutsCoincide(t *testing.T) {
	// DpnII (^GATC) and MboI (^GATC) cut at SAME position → two 0-length frags, no bridging.
	eA := enzyme.DB["DpnII"]
	eB := enzyme.DB["MboI"]
	seq := []byte("AAAGATCAAAGATC")
	frags := Digest(seq, []enzyme.Enzyme{eA, eB}, 0, 1<<30)
	if len(frags) != 2 {
		t.Fatalf("expect 2 zero-length frags (one per site), got %d (%v)", len(frags), frags)
	}
	if frags[0].End-frags[0].Start != 0 || frags[1].End-frags[1].Start != 0 {
		t.Fatalf("expected zero-length fragments, got %v", frags)
	}
}

func TestCutOffsetFallback_NoCaretUsesMidsite(t *testing.T) {
	// Non-strict plan keeps mid-site fallback behavior.
	fake := enzyme.Enzyme{Name: "Fake", Recognition: "AAAA", CutIndex: 0}
	frags := Digest([]byte("AAAAAA"), []enzyme.Enzyme{fake}, 1, 1<<30)
	if len(frags) == 0 {
		t.Fatalf("expected at least one fragment with mid-site fallback")
	}
}
