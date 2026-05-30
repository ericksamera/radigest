package digest

import (
	"testing"

	"github.com/ericksamera/radigest/internal/enzyme"
)

// toy chr:---GAATTC---TTAA---GAATTC---
// EcoRI cuts at index 1 within GAATTC (after 'G'): positions 5 and 16 in toyChr.
const toyChr = "AAAAGAATTCTTAAAGAATTC"

func TestSingleDigest_ConsecutiveCuts(t *testing.T) {
	eA := enzyme.DB["EcoRI"]
	frags := Digest([]byte(toyChr), []enzyme.Enzyme{eA}, 0, 1<<30)
	if len(frags) != 1 {
		t.Fatalf("want 1 fragment, got %d: %#v", len(frags), frags)
	}
	want := Fragment{Start: 5, End: 16}
	if frags[0] != want {
		t.Fatalf("mismatch: got %+v want %+v", frags[0], want)
	}
}

func TestDoubleDigest_AB_BA_Filter(t *testing.T) {
	eA := enzyme.DB["EcoRI"]
	eB := enzyme.DB["MseI"]

	frags := Digest([]byte(toyChr), []enzyme.Enzyme{eA, eB}, 0, 1<<30)
	if len(frags) != 2 { // should keep A-B and B-A only
		t.Fatalf("want 2 kept fragments, got %d: %#v", len(frags), frags)
	}
	want := []Fragment{
		{Start: 5, End: 11},  // EcoRI-MseI
		{Start: 11, End: 16}, // MseI-EcoRI
	}
	for i := range want {
		if frags[i] != want[i] {
			t.Fatalf("frag %d mismatch: got %+v want %+v", i, frags[i], want[i])
		}
	}
}

func TestExactScannerDetectsOverlappingMotifs(t *testing.T) {
	eA := enzyme.Enzyme{Name: "FakeExact", Recognition: "AA^A"}
	seq := []byte("AAAAA")

	frags := Digest(seq, []enzyme.Enzyme{eA}, 0, 1<<30)
	want := []Fragment{
		{Start: 2, End: 3},
		{Start: 3, End: 4},
	}
	if len(frags) != len(want) {
		t.Fatalf("got %d fragments, want %d: %#v", len(frags), len(want), frags)
	}
	for i := range want {
		if frags[i] != want[i] {
			t.Fatalf("fragment %d got %+v, want %+v", i, frags[i], want[i])
		}
	}
}

func TestDigestStatsMatchesDigest(t *testing.T) {
	cases := []struct {
		name string
		seq  []byte
		ens  []enzyme.Enzyme
		opt  Options
		min  int
		max  int
	}{
		{
			name: "single internal",
			seq:  []byte("AATTCCCCAATT"),
			ens:  []enzyme.Enzyme{{Name: "MluCI", Recognition: "^AATT"}},
			min:  1,
			max:  100,
		},
		{
			name: "double AB BA",
			seq:  []byte("AAAAGAATTCTTAAAGAATTC"),
			ens: []enzyme.Enzyme{
				{Name: "EcoRI", Recognition: "G^AATTC"},
				{Name: "MseI", Recognition: "T^TAA"},
			},
			min: 1,
			max: 100,
		},
		{
			name: "double include ends",
			seq:  []byte("AAAAGAATTCTTAAAGAATTC"),
			ens: []enzyme.Enzyme{
				{Name: "EcoRI", Recognition: "G^AATTC"},
				{Name: "MseI", Recognition: "T^TAA"},
			},
			opt: Options{IncludeEnds: true},
			min: 1,
			max: 100,
		},
		{
			name: "same enzyme allowed",
			seq:  []byte("AAAAGAATTCAAAAGAATTCAAA"),
			ens: []enzyme.Enzyme{
				{Name: "EcoRI", Recognition: "G^AATTC"},
				{Name: "NcoI", Recognition: "C^CATGG"},
			},
			opt: Options{AllowSame: true},
			min: 1,
			max: 100,
		},
		{
			name: "coincident zero length",
			seq:  []byte("AAAGATCAAAGATC"),
			ens: []enzyme.Enzyme{
				{Name: "DpnII", Recognition: "^GATC"},
				{Name: "MboI", Recognition: "^GATC"},
			},
			min: 0,
			max: 100,
		},
		{
			name: "degenerate motif",
			seq:  []byte("AAGCAGCAAGCTGC"),
			ens:  []enzyme.Enzyme{{Name: "ApeKI", Recognition: "G^CWGC"}},
			min:  1,
			max:  100,
		},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			plan := NewPlanWithOptions(tc.ens, tc.opt)
			frags := plan.Digest(tc.seq, tc.min, tc.max)
			stats := plan.DigestStats(tc.seq, tc.min, tc.max)

			var wantBases int
			for _, fr := range frags {
				wantBases += fr.End - fr.Start
			}
			if stats.Fragments != len(frags) || stats.Bases != wantBases {
				t.Fatalf(
					"DigestStats=(frags=%d bases=%d), Digest=(frags=%d bases=%d)",
					stats.Fragments,
					stats.Bases,
					len(frags),
					wantBases,
				)
			}
		})
	}
}
