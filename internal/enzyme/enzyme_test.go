package enzyme

import "testing"

func TestLookup(t *testing.T) {
	en := DB["EcoRI"]
	if en.Name != "EcoRI" || en.Recognition != "G^AATTC" {
		t.Fatalf("EcoRI not found or wrong: %+v", en)
	}
	if _, ok := DB["ImaginaryI"]; ok {
		t.Fatalf("unexpected fake enzyme present")
	}
}
