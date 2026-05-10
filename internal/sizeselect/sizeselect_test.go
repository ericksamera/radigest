package sizeselect

import (
	"math"
	"testing"
)

func TestHardWeight(t *testing.T) {
	sel, err := New(Config{Model: ModelHard, Min: 100, Max: 200, ScoreMin: 1, ScoreMax: 500})
	if err != nil {
		t.Fatal(err)
	}
	if sel.Weight(99) != 0 || sel.Weight(100) != 1 || sel.Weight(200) != 1 || sel.Weight(201) != 0 {
		t.Fatalf("hard weights not as expected")
	}
}

func TestNormalWeight(t *testing.T) {
	sel, err := New(Config{Model: ModelNormal, Min: 100, Max: 200, ScoreMin: 1, ScoreMax: 500, Mean: 150, SD: 25})
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(sel.Weight(150)-1) > 1e-12 {
		t.Fatalf("normal at mean should be 1, got %g", sel.Weight(150))
	}
	if math.Abs(sel.Weight(125)-sel.Weight(175)) > 1e-12 {
		t.Fatalf("normal weights should be symmetric")
	}
}

func TestNormalWeightAllowsMeanOutsideHardWindow(t *testing.T) {
	sel, err := New(Config{Model: ModelNormal, Min: 300, Max: 600, ScoreMin: 1, ScoreMax: 2000, Mean: 275, SD: 85})
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(sel.Weight(275)-1) > 1e-12 {
		t.Fatalf("normal at empirical mean should be 1, got %g", sel.Weight(275))
	}
	if !sel.InScoreRange(50) || sel.InHardWindow(50) {
		t.Fatalf("score/hard windows not handled as expected")
	}
}

func TestTriangularWeight(t *testing.T) {
	sel, err := New(Config{Model: ModelTriangular, Min: 100, Max: 200, ScoreMin: 1, ScoreMax: 500, Mean: 150})
	if err != nil {
		t.Fatal(err)
	}
	if sel.Weight(100) != 0 || sel.Weight(150) != 1 || sel.Weight(200) != 0 {
		t.Fatalf("triangular weights unexpected: 100=%g 150=%g 200=%g", sel.Weight(100), sel.Weight(150), sel.Weight(200))
	}
	if math.Abs(sel.Weight(125)-0.5) > 1e-12 || math.Abs(sel.Weight(175)-0.5) > 1e-12 {
		t.Fatalf("triangular shoulder weights unexpected")
	}
}

func TestSoftWindowWeight(t *testing.T) {
	sel, err := New(Config{Model: ModelSoftWindow, Min: 100, Max: 200, ScoreMin: 1, ScoreMax: 500, EdgeSD: 10})
	if err != nil {
		t.Fatal(err)
	}
	if !(sel.Weight(150) > sel.Weight(100) && sel.Weight(100) > sel.Weight(50)) {
		t.Fatalf("soft-window weights should rise into window")
	}
	if !(sel.Weight(150) > sel.Weight(200) && sel.Weight(200) > sel.Weight(250)) {
		t.Fatalf("soft-window weights should fall out of window")
	}
}

func TestInvalidConfigs(t *testing.T) {
	bad := []Config{
		{Model: "nope", Min: 1, Max: 10, ScoreMin: 1, ScoreMax: 10},
		{Model: ModelNormal, Min: 1, Max: 10, ScoreMin: 1, ScoreMax: 10, Mean: 5, SD: 0},
		{Model: ModelSoftWindow, Min: 1, Max: 10, ScoreMin: 1, ScoreMax: 10, EdgeSD: -1},
		{Model: ModelTriangular, Min: 1, Max: 10, ScoreMin: 1, ScoreMax: 10, Mean: 10},
		{Model: ModelHard, Min: 1, Max: 10, ScoreMin: 20, ScoreMax: 10},
	}
	for _, cfg := range bad {
		if _, err := New(cfg); err == nil {
			t.Fatalf("New(%+v) returned nil error", cfg)
		}
	}
}

func TestStatsAdd(t *testing.T) {
	sel, err := New(Config{Model: ModelHard, Min: 100, Max: 200, ScoreMin: 1, ScoreMax: 500})
	if err != nil {
		t.Fatal(err)
	}
	st := NewStats(sel)
	st.Add(120, true, 1)
	st.Add(80, false, 0)
	if st.RawFragmentsScored != 2 || st.RawBasesScored != 200 || st.RawFragmentsInWindow != 1 || st.RawBasesInWindow != 120 {
		t.Fatalf("raw stats wrong: %+v", st)
	}
	if st.WeightedFragments != 1 || st.WeightedBases != 120 || st.MeanWeightedLength != 120 {
		t.Fatalf("weighted stats wrong: %+v", st)
	}
}
