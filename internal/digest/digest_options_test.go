package digest

import (
	"reflect"
	"testing"

	"github.com/ericksamera/radigest/internal/enzyme"
)

func TestAllowSame_EnablesAAFragments(t *testing.T) {
	eA := enzyme.DB["EcoRI"]
	eB := enzyme.DB["NcoI"] // absent
	seq := []byte("AAAAGAATTCAAAAGAATTCAAA")

	pNo := NewPlanWithOptions([]enzyme.Enzyme{eA, eB}, Options{AllowSame: false})
	if got := pNo.Digest(seq, 1, 1<<30); len(got) != 0 {
		t.Fatalf("default should drop AA/BB, got %d", len(got))
	}
	pYes := NewPlanWithOptions([]enzyme.Enzyme{eA, eB}, Options{AllowSame: true})
	if got := pYes.Digest(seq, 1, 1<<30); len(got) == 0 {
		t.Fatalf("AllowSame should keep AA/BB, got 0")
	}
}

func TestIncludeEnds_EnablesTerminalFragmentsSingleDigest(t *testing.T) {
	eA := enzyme.DB["MluCI"] // ^AATT
	seq := []byte("CCCCAATTGGGGAATTTT")

	pNo := NewPlanWithOptions([]enzyme.Enzyme{eA}, Options{IncludeEnds: false})
	if got, want := pNo.Digest(seq, 1, 1<<30), ([]Fragment{{Start: 4, End: 12}}); !reflect.DeepEqual(got, want) {
		t.Fatalf("default fragments mismatch: got %#v want %#v", got, want)
	}

	pYes := NewPlanWithOptions([]enzyme.Enzyme{eA}, Options{IncludeEnds: true})
	got := pYes.Digest(seq, 1, 1<<30)
	want := []Fragment{{Start: 0, End: 4}, {Start: 4, End: 12}, {Start: 12, End: 18}}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("include-ends fragments mismatch: got %#v want %#v", got, want)
	}
}

func TestIncludeEnds_EnablesTerminalFragmentsDoubleDigest(t *testing.T) {
	eA := enzyme.DB["EcoRI"]
	eB := enzyme.DB["MseI"]
	p := NewPlanWithOptions([]enzyme.Enzyme{eA, eB}, Options{IncludeEnds: true})

	got := p.Digest([]byte(toyChr), 1, 1<<30)
	want := []Fragment{
		{Start: 0, End: 5},
		{Start: 5, End: 11},
		{Start: 11, End: 16},
		{Start: 16, End: len(toyChr)},
	}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("include-ends fragments mismatch: got %#v want %#v", got, want)
	}
}

func TestStrictCuts_PanicsOnFallback(t *testing.T) {
	defer func() {
		if r := recover(); r == nil {
			t.Fatalf("StrictCuts should panic when fallback would be used")
		}
	}()
	fake := enzyme.Enzyme{Name: "Fake", Recognition: "AAAA", CutIndex: 0} // no caret, CutIndex==0
	_ = NewPlanWithOptions([]enzyme.Enzyme{fake}, Options{StrictCuts: true})
}
