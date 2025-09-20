package sim

import (
	"bytes"
	"math"
	"testing"
)

func gcFrac(b []byte) float64 {
	if len(b) == 0 { return 0 }
	gc := 0
	for _, x := range b {
		if x == 'G' || x == 'C' { gc++ }
	}
	return float64(gc) / float64(len(b))
}

func TestMake_LengthAndGC(t *testing.T) {
	N := 10000
	seq := Make(N, 0.42, 123)
	if len(seq) != N {
		t.Fatalf("length: got %d want %d", len(seq), N)
	}
	got := gcFrac(seq)
	want := 0.42
	tol := 0.5/float64(N) + 1e-12 // nearest-integer rounding
	if math.Abs(got-want) > tol {
		t.Fatalf("gc: got %.6f want %.6f (tol %.6g)", got, want, tol)
	}
}

func TestMake_SeedDeterministic(t *testing.T) {
	a := Make(5000, 0.50, 42)
	b := Make(5000, 0.50, 42)
	if !bytes.Equal(a, b) {
		t.Fatalf("same seed should reproduce sequence")
	}
	c := Make(5000, 0.50, 43)
	if bytes.Equal(a, c) {
		t.Fatalf("different seed unexpectedly produced identical sequence")
	}
}

func TestMake_GCExtremesAndClamp(t *testing.T) {
	at := Make(1000, 0, 7)
	for _, x := range at {
		if x == 'G' || x == 'C' {
			t.Fatalf("expected only A/T when gc=0, saw %c", x)
		}
	}
	gc := Make(1000, 1, 7)
	for _, x := range gc {
		if x == 'A' || x == 'T' {
			t.Fatalf("expected only G/C when gc=1, saw %c", x)
		}
	}
	if len(Make(0, 0.5, 1)) != 0 {
		t.Fatalf("length zero should return empty slice")
	}
	// clamp checks
	if g := gcFrac(Make(100, -0.1, 1)); g != 0 {
		t.Fatalf("gc clamp low failed: got %.3f", g)
	}
	if g := gcFrac(Make(100, 1.5, 1)); g != 1 {
		t.Fatalf("gc clamp high failed: got %.3f", g)
	}
}
