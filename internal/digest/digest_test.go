package digest

import (
	"testing"

	"radigest/internal/enzyme"
)

// toy chr:---GAATTC---TTAA---GAATTC---
// EcoRI cuts 6→^; MseI cuts 5→^
const toyChr = "AAAAGAATTCTTAAAGAATTC"

func TestDoubleDigest_AB_BA_Filter(t *testing.T) {
	eA := enzyme.DB["EcoRI"]
	eB := enzyme.DB["MseI"]

	frags := Digest([]byte(toyChr), []enzyme.Enzyme{eA, eB}, 0, 1<<30)
	if len(frags) != 2 { // should keep A-B and B-A only
		t.Fatalf("want 2 kept fragments, got %d: %#v", len(frags), frags)
	}
	// check coordinates (1-based end-exclusive for this test)
	want := []Fragment{
		{Start: 5, End: 11}, // EcoRI-MseI
		{Start: 11, End: 16}, // MseI-EcoRI
	}
	for i := range want {
		if frags[i] != want[i] {
			t.Fatalf("frag %d mismatch: got %+v want %+v", i, frags[i], want[i])
		}
	}
}
