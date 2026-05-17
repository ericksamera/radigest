package fragmentfasta

import (
	"os"
	"strings"
	"testing"

	"github.com/ericksamera/radigest/internal/digest"
)

func TestWriter(t *testing.T) {
	tmp, err := os.CreateTemp("", "fragments*.fa")
	if err != nil {
		t.Fatal(err)
	}
	path := tmp.Name()
	_ = tmp.Close()
	defer func() { _ = os.Remove(path) }()

	w, err := New(path)
	if err != nil {
		t.Fatal(err)
	}
	seq := []byte(strings.Repeat("ACGT", 25))
	if err := w.Write("chr 1", 1, digest.Fragment{Start: 0, End: 85}, seq); err != nil {
		t.Fatal(err)
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}

	raw, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	text := string(raw)
	if !strings.HasPrefix(text, ">chr%201_1 chrom=chr%201 start0=0 end0=85 length=85\n") {
		t.Fatalf("unexpected header/body: %q", text)
	}
	lines := strings.Split(strings.TrimSpace(text), "\n")
	if len(lines) != 3 {
		t.Fatalf("expected header + 2 wrapped sequence lines, got %d: %q", len(lines), text)
	}
	if len(lines[1]) != 80 || len(lines[2]) != 5 {
		t.Fatalf("unexpected sequence wrapping: line lengths %d and %d", len(lines[1]), len(lines[2]))
	}
}

func TestWriterRejectsInvalidFragment(t *testing.T) {
	tmp, err := os.CreateTemp("", "fragments*.fa")
	if err != nil {
		t.Fatal(err)
	}
	path := tmp.Name()
	_ = tmp.Close()
	defer func() { _ = os.Remove(path) }()

	w, err := New(path)
	if err != nil {
		t.Fatal(err)
	}
	defer func() { _ = w.Close() }()

	if err := w.Write("chr1", 1, digest.Fragment{Start: 0, End: 10}, []byte("ACGT")); err == nil {
		t.Fatalf("expected invalid-fragment error")
	}
}

func TestDisabledWriterNoops(t *testing.T) {
	w, err := New("")
	if err != nil {
		t.Fatal(err)
	}
	if err := w.Write("chr1", 1, digest.Fragment{Start: 0, End: 1}, []byte("A")); err != nil {
		t.Fatal(err)
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
}

func TestNilWriterNoops(t *testing.T) {
	var w *Writer
	if err := w.Write("chr1", 1, digest.Fragment{Start: 0, End: 1}, []byte("A")); err != nil {
		t.Fatal(err)
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
}
