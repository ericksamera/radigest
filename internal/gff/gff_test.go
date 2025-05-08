package gff

import (
	"bytes"
	"strings"
	"testing"

	"radigest/internal/digest"
)

func TestWrite(t *testing.T) {
	buf := &bytes.Buffer{}
	frags := []digest.Fragment{{Start: 4, End: 10}} // 0-based half-open
	if err := Write(buf, "chr1", frags); err != nil {
		t.Fatal(err)
	}
	want := "##gff-version 3\nchr1\tradigest\tfragment\t5\t10\t.\t+\t.\tID=frag1\n"
	got := buf.String()
	if strings.TrimSpace(got) != strings.TrimSpace(want) {
		t.Fatalf("mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}
