package digest

import (
	"testing"

	"radigest/internal/enzyme"
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

func TestStrictCuts_PanicsOnFallback(t *testing.T) {
	defer func() {
		if r := recover(); r == nil {
			t.Fatalf("StrictCuts should panic when fallback would be used")
		}
	}()
	fake := enzyme.Enzyme{Name: "Fake", Recognition: "AAAA", CutIndex: 0} // no caret, CutIndex==0
	_ = NewPlanWithOptions([]enzyme.Enzyme{fake}, Options{StrictCuts: true})
}
