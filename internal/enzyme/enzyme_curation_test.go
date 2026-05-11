package enzyme

import (
	"encoding/json"
	"os"
	"testing"
)

type unsupportedEnzyme struct {
	Name string `json:"name"`
	Site string `json:"site"`
}

type aliasGroup struct {
	Canonical string   `json:"canonical"`
	Aliases   []string `json:"aliases"`
	Reason    string   `json:"reason"`
}

func TestUnsupportedEnzymesAreDocumentedAndAbsentFromDefaultDB(t *testing.T) {
	raw, err := os.ReadFile("enzymes.unsupported.json")
	if err != nil {
		t.Fatalf("read unsupported enzyme metadata: %v", err)
	}

	var unsupported []unsupportedEnzyme
	if err := json.Unmarshal(raw, &unsupported); err != nil {
		t.Fatalf("parse unsupported enzyme metadata: %v", err)
	}

	if len(unsupported) == 0 {
		t.Fatalf("unsupported enzyme metadata is empty")
	}

	seen := make(map[string]bool)
	for _, e := range unsupported {
		if e.Name == "" {
			t.Fatalf("unsupported enzyme has empty name: %+v", e)
		}
		if e.Site == "" {
			t.Fatalf("unsupported enzyme %s has empty site", e.Name)
		}
		if _, ok := DB[e.Name]; ok {
			t.Fatalf("unsupported enzyme %s is still present in default DB", e.Name)
		}
		seen[e.Name] = true
	}

	expected := []string{
		"AciI",
		"BbvCI",
		"Bpu10I",
		"BseYI",
		"BsmAI",
		"BsmI",
		"BssSI-v2",
		"PleI",
		"I-CeuI",
		"I-SceI",
		"PI-PspI",
		"PI-SceI",
		"SfiI",
		"SgrAI",
		"NgoMIV",
	}
	for _, name := range expected {
		if !seen[name] {
			t.Fatalf("expected unsupported enzyme %s is not documented", name)
		}
	}
}

func TestAliasMetadataReferencesSupportedEquivalentEnzymes(t *testing.T) {
	raw, err := os.ReadFile("enzymes.aliases.json")
	if err != nil {
		t.Fatalf("read enzyme alias metadata: %v", err)
	}

	var groups []aliasGroup
	if err := json.Unmarshal(raw, &groups); err != nil {
		t.Fatalf("parse enzyme alias metadata: %v", err)
	}

	if len(groups) == 0 {
		t.Fatalf("enzyme alias metadata is empty")
	}

	for _, group := range groups {
		canonical, ok := DB[group.Canonical]
		if !ok {
			t.Fatalf("canonical alias enzyme %s is not present in default DB", group.Canonical)
		}
		if len(group.Aliases) == 0 {
			t.Fatalf("alias group %s has no aliases", group.Canonical)
		}

		for _, aliasName := range group.Aliases {
			alias, ok := DB[aliasName]
			if !ok {
				t.Fatalf("alias enzyme %s is not present in default DB", aliasName)
			}
			if alias.Recognition != canonical.Recognition || alias.CutIndex != canonical.CutIndex {
				t.Fatalf(
					"alias %s is not sequence-equivalent to %s: alias=%+v canonical=%+v",
					aliasName,
					group.Canonical,
					alias,
					canonical,
				)
			}
		}
	}
}
