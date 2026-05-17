package digest

import (
	"reflect"
	"strings"
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

func TestTryNewPlanWithOptions_StrictCutsReturnsErrorOnFallback(t *testing.T) {
	fake := enzyme.Enzyme{Name: "Fake", Recognition: "AAAA", CutIndex: 0} // no caret, CutIndex==0
	if _, err := TryNewPlanWithOptions([]enzyme.Enzyme{fake}, Options{StrictCuts: true}); err == nil {
		t.Fatalf("StrictCuts should return an error when fallback would be used")
	}
}

func TestTryNewPlanWithOptions_RejectsInvalidIUPAC(t *testing.T) {
	fake := enzyme.Enzyme{Name: "Bad", Recognition: "A^X", CutIndex: 1}
	_, err := TryNewPlanWithOptions([]enzyme.Enzyme{fake}, Options{})
	if err == nil {
		t.Fatalf("TryNewPlanWithOptions returned nil error for invalid IUPAC symbol")
	}
	if !strings.Contains(err.Error(), "invalid IUPAC") {
		t.Fatalf("expected invalid IUPAC error, got %v", err)
	}
}

func TestTryNewPlanWithOptions_ValidEnzymesCompile(t *testing.T) {
	plan, err := TryNewPlanWithOptions([]enzyme.Enzyme{enzyme.DB["EcoRI"], enzyme.DB["MseI"]}, Options{})
	if err != nil {
		t.Fatalf("TryNewPlanWithOptions returned error for valid enzymes: %v", err)
	}
	if got := plan.Digest([]byte(toyChr), 1, 1<<30); len(got) == 0 {
		t.Fatalf("valid plan produced no fragments")
	}
}
