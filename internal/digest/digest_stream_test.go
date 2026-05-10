package digest

import (
	"errors"
	"reflect"
	"testing"

	"github.com/KPU-AGC/radigest/internal/enzyme"
)

func collectDigestEach(t *testing.T, p Plan, seq []byte, min, max int) []Fragment {
	t.Helper()
	var got []Fragment
	if err := p.DigestEach(seq, min, max, func(fr Fragment) error {
		got = append(got, fr)
		return nil
	}); err != nil {
		t.Fatalf("DigestEach returned error: %v", err)
	}
	return got
}

func TestDigestEachMatchesExpectedSingle(t *testing.T) {
	p := NewPlan([]enzyme.Enzyme{enzyme.DB["EcoRI"]})
	got := collectDigestEach(t, p, []byte(toyChr), 0, 1<<30)
	want := []Fragment{{Start: 5, End: 16}}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("fragments mismatch: got %#v want %#v", got, want)
	}
}

func TestDigestEachMatchesExpectedDouble(t *testing.T) {
	p := NewPlan([]enzyme.Enzyme{enzyme.DB["EcoRI"], enzyme.DB["MseI"]})
	got := collectDigestEach(t, p, []byte(toyChr), 0, 1<<30)
	want := []Fragment{{Start: 5, End: 11}, {Start: 11, End: 16}}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("fragments mismatch: got %#v want %#v", got, want)
	}
}

func TestDigestEachDoubleDigestOrdersByCutCoordinate(t *testing.T) {
	// A is recognized first by window start but cuts later than B. This catches
	// implementations that scan by recognition-position order instead of cut order.
	eA := enzyme.Enzyme{Name: "A", Recognition: "AACC^"} // window start 0, cut 4
	eB := enzyme.Enzyme{Name: "B", Recognition: "^CC"}   // window start 2, cut 2
	p := NewPlan([]enzyme.Enzyme{eA, eB})

	got := collectDigestEach(t, p, []byte("AACC"), 0, 1<<30)
	want := []Fragment{{Start: 2, End: 4}}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("fragments mismatch: got %#v want %#v", got, want)
	}
}

func TestDigestEachPropagatesEmitError(t *testing.T) {
	p := NewPlan([]enzyme.Enzyme{enzyme.DB["EcoRI"], enzyme.DB["MseI"]})
	wantErr := errors.New("stop")
	calls := 0
	err := p.DigestEach([]byte(toyChr), 0, 1<<30, func(Fragment) error {
		calls++
		return wantErr
	})
	if !errors.Is(err, wantErr) {
		t.Fatalf("got error %v want %v", err, wantErr)
	}
	if calls != 1 {
		t.Fatalf("emit called %d times, want 1", calls)
	}
}

func TestDigestEachRejectsNilEmit(t *testing.T) {
	p := NewPlan([]enzyme.Enzyme{enzyme.DB["EcoRI"]})
	if err := p.DigestEach([]byte(toyChr), 0, 1<<30, nil); err == nil {
		t.Fatal("DigestEach with nil emit returned nil error")
	}
}
