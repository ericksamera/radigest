package enzyme

import "testing"

func TestIUPAC_N_MatchesAny(t *testing.T) {
	p, err := CompilePattern("N")
	if err != nil { t.Fatal(err) }
	for _, b := range []byte("ACGT") {
		if !Match(p, []byte{b}) {
			t.Fatalf("N should match %c", b)
		}
	}
}
