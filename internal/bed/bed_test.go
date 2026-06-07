package bed

import (
	"bytes"
	"os"
	"strings"
	"testing"

	"github.com/ericksamera/radigest/internal/digest"
)

func TestWriter(t *testing.T) {
	tmp, err := os.CreateTemp("", "fragments*.bed")
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
	if err := w.Write("chr1", 1, digest.Fragment{Start: 10, End: 25}); err != nil {
		t.Fatal(err)
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}

	raw, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	want := "chr1\t10\t25\tchr1_1\t0\t+\n"
	if string(raw) != want {
		t.Fatalf("BED mismatch\nwant: %q\ngot:  %q", want, string(raw))
	}
}

func TestWriterEscapesFragmentName(t *testing.T) {
	tmp, err := os.CreateTemp("", "fragments*.bed")
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
	if err := w.Write("chr|1%bad", 2, digest.Fragment{Start: 0, End: 5}); err != nil {
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
	if !strings.Contains(text, "\tchr|1%25bad_2\t") {
		t.Fatalf("fragment name was not escaped: %q", text)
	}
}

func TestDisabledWriterNoops(t *testing.T) {
	w, err := New("")
	if err != nil {
		t.Fatal(err)
	}
	if err := w.Write("chr1", 1, digest.Fragment{Start: 0, End: 1}); err != nil {
		t.Fatal(err)
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
}

func TestNilWriterNoops(t *testing.T) {
	var w *Writer
	if err := w.Write("chr1", 1, digest.Fragment{Start: 0, End: 1}); err != nil {
		t.Fatal(err)
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
}

func TestWriterRejectsInvalidFragment(t *testing.T) {
	var buf bytes.Buffer
	w, err := NewTo("-", &buf)
	if err != nil {
		t.Fatal(err)
	}
	defer func() { _ = w.Close() }()
	if err := w.Write("chr1", 1, digest.Fragment{Start: 10, End: 5}); err == nil {
		t.Fatalf("expected invalid-fragment error")
	}
}
