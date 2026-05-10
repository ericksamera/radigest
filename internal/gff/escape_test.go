package gff

import "testing"

func TestEscapeSeqID(t *testing.T) {
	got := EscapeSeqID("chr 1;bad=2,50%\t")
	want := "chr%201%3Bbad%3D2%2C50%25%09"
	if got != want {
		t.Fatalf("EscapeSeqID mismatch: got %q want %q", got, want)
	}
}

func TestEscapeAttributeValue(t *testing.T) {
	got := EscapeAttributeValue("chr 1;bad=2,50%\t")
	want := "chr%201%3Bbad%3D2%2C50%25%09"
	if got != want {
		t.Fatalf("EscapeAttributeValue mismatch: got %q want %q", got, want)
	}
}

func TestFragmentAttributesEscapesID(t *testing.T) {
	got := FragmentAttributes("chr 1;bad=2,50%\t", 7, 42)
	want := "ID=chr%201%3Bbad%3D2%2C50%25%09_7;Length=42"
	if got != want {
		t.Fatalf("FragmentAttributes mismatch: got %q want %q", got, want)
	}
}
