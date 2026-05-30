package digest

import (
	"errors"
	"reflect"
	"testing"

	"github.com/ericksamera/radigest/internal/enzyme"
)

func fragmentsEqual(a, b []Fragment) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

func collectCutsForEnzyme(t *testing.T, e enzyme.Enzyme, seq []byte) []int {
	t.Helper()
	plan := NewPlan([]enzyme.Enzyme{e})
	return plan.Cuts(seq)
}

func collectDigestCutsEach(t *testing.T, cutsA, cutsB []int, seqLen, min, max int, opt Options) []Fragment {
	t.Helper()
	var got []Fragment
	if err := DigestCutsEach(cutsA, cutsB, seqLen, min, max, opt, func(fr Fragment) error {
		got = append(got, fr)
		return nil
	}); err != nil {
		t.Fatalf("DigestCutsEach returned error: %v", err)
	}
	return got
}

func TestCutsReturnsSortedCutCoordinates(t *testing.T) {
	got := collectCutsForEnzyme(t, enzyme.DB["EcoRI"], []byte(toyChr))
	want := []int{5, 16}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("EcoRI cuts mismatch: got %#v want %#v", got, want)
	}
}

func TestCutsSupportsDegenerateMotifsAndReferenceNPolicy(t *testing.T) {
	apekiCuts := collectCutsForEnzyme(t, enzyme.Enzyme{Name: "ApeKI", Recognition: "G^CWGC"}, []byte("AAGCAGCAAGCTGC"))
	if want := []int{3, 10}; !reflect.DeepEqual(apekiCuts, want) {
		t.Fatalf("ApeKI cuts mismatch: got %#v want %#v", apekiCuts, want)
	}

	hinfCuts := collectCutsForEnzyme(t, enzyme.Enzyme{Name: "HinfI", Recognition: "G^ANTC"}, []byte("GAATCCCCGANTCCCCGAATC"))
	if want := []int{1, 17}; !reflect.DeepEqual(hinfCuts, want) {
		t.Fatalf("HinfI cuts mismatch: got %#v want %#v", hinfCuts, want)
	}
}

func TestCutsEachPropagatesEmitError(t *testing.T) {
	wantErr := errors.New("stop")
	plan := NewPlan([]enzyme.Enzyme{enzyme.DB["EcoRI"]})
	err := plan.CutsEach([]byte(toyChr), func(int) error { return wantErr })
	if !errors.Is(err, wantErr) {
		t.Fatalf("CutsEach error = %v, want %v", err, wantErr)
	}
}

func TestCutsEachRejectsNilEmit(t *testing.T) {
	plan := NewPlan([]enzyme.Enzyme{enzyme.DB["EcoRI"]})
	if err := plan.CutsEach([]byte(toyChr), nil); err == nil {
		t.Fatal("CutsEach with nil emit returned nil error")
	}
}

func TestDigestCutsEachRejectsNilEmit(t *testing.T) {
	if err := DigestCutsEach([]int{1, 5}, nil, 10, 1, 100, Options{}, nil); err == nil {
		t.Fatal("DigestCutsEach with nil emit returned nil error")
	}
}

func TestDigestCutsEachMatchesPlanDigestEach(t *testing.T) {
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
			name: "single include ends no cut",
			seq:  []byte("CCCCCCCC"),
			ens:  []enzyme.Enzyme{{Name: "EcoRI", Recognition: "G^AATTC"}},
			opt:  Options{IncludeEnds: true},
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
			name: "same enzyme suppressed",
			seq:  []byte("AAAAGAATTCAAAAGAATTCAAA"),
			ens: []enzyme.Enzyme{
				{Name: "EcoRI", Recognition: "G^AATTC"},
				{Name: "NcoI", Recognition: "C^CATGG"},
			},
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
			name: "double second enzyme has no cuts",
			seq:  []byte("AAAAGAATTCAAAAGAATTCAAA"),
			ens: []enzyme.Enzyme{
				{Name: "EcoRI", Recognition: "G^AATTC"},
				{Name: "NcoI", Recognition: "C^CATGG"},
			},
			opt: Options{IncludeEnds: true},
			min: 1,
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
			want := collectDigestEach(t, plan, tc.seq, tc.min, tc.max)

			cutsA := collectCutsForEnzyme(t, tc.ens[0], tc.seq)
			var cutsB []int
			if len(tc.ens) > 1 {
				cutsB = collectCutsForEnzyme(t, tc.ens[1], tc.seq)
			}

			got := collectDigestCutsEach(t, cutsA, cutsB, len(tc.seq), tc.min, tc.max, tc.opt)
			if !fragmentsEqual(got, want) {
				t.Fatalf("DigestCutsEach mismatch:\n got %#v\nwant %#v", got, want)
			}

			gotSlice := DigestCuts(cutsA, cutsB, len(tc.seq), tc.min, tc.max, tc.opt)
			if !fragmentsEqual(gotSlice, want) {
				t.Fatalf("DigestCuts mismatch:\n got %#v\nwant %#v", gotSlice, want)
			}
		})
	}
}
