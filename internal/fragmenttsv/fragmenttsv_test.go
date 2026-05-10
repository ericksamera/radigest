package fragmenttsv

import (
	"os"
	"strings"
	"testing"

	"github.com/KPU-AGC/radigest/internal/digest"
)

func TestWriter(t *testing.T) {
	tmp, err := os.CreateTemp("", "fragments*.tsv")
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
	if err := w.Write("chr1", digest.Fragment{Start: 10, End: 25}, true, 0.75); err != nil {
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
	if !strings.HasPrefix(text, "chrom\tstart0\tend0\tlength\thard_kept\tsize_weight\n") {
		t.Fatalf("missing header: %q", text)
	}
	if !strings.Contains(text, "chr1\t10\t25\t15\ttrue\t0.75\n") {
		t.Fatalf("unexpected body: %q", text)
	}
}

func TestDisabledWriterNoops(t *testing.T) {
	w, err := New("")
	if err != nil {
		t.Fatal(err)
	}
	if err := w.Write("chr1", digest.Fragment{Start: 0, End: 1}, false, 0); err != nil {
		t.Fatal(err)
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
}

func TestNilWriterNoops(t *testing.T) {
	var w *Writer
	if err := w.Write("chr1", digest.Fragment{Start: 0, End: 1}, false, 0); err != nil {
		t.Fatal(err)
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
}
