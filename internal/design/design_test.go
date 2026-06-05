package design

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/ericksamera/radigest/internal/screen"
	"github.com/ericksamera/radigest/internal/sizeselect"
)

func TestEvaluateSummaryFeasibleCandidate(t *testing.T) {
	summary := screen.PairSummary{
		Enzymes: []string{"EcoRI", "MseI"},
		SizeSelection: sizeselect.Stats{
			WeightedBases:        2500,
			WeightedFragments:    100,
			MeanWeightedLength:   400,
			RawBasesInWindow:     2500,
			RawFragmentsInWindow: 100,
		},
		Screening: screen.ScreeningStats{Records: 1, CachedCutSites: 10, CacheMemoryEstimateBytes: 80},
	}
	budget := SequencingBudget{ReadLayout: "pe", ReadLength: 150, LaneReadPairs: 2000, Lanes: 1, UsableReadFraction: 1, Samples: 1, DesiredDepth: 10}
	target := DesignTarget{TargetGenomePct: 2.5, CoverageTolerancePct: 0.01, Objective: ObjectiveBalanced}

	cand := EvaluateSummary(summary, 100000, budget, target, DefaultScoreWeights())
	if !cand.Feasible {
		t.Fatalf("candidate should be feasible: %+v", cand)
	}
	if cand.GeneratedWeightedGenomePct != 2.5 {
		t.Fatalf("weighted genome pct = %g, want 2.5", cand.GeneratedWeightedGenomePct)
	}
	if cand.ExpectedMeanDepth != 20 {
		t.Fatalf("expected depth = %g, want 20", cand.ExpectedMeanDepth)
	}
	if cand.MaxSamplesTotalFullTarget != 2 {
		t.Fatalf("max samples total = %d, want 2", cand.MaxSamplesTotalFullTarget)
	}
}

func TestEvaluateSummaryDepthShortfall(t *testing.T) {
	summary := screen.PairSummary{
		Enzymes:       []string{"EcoRI", "MseI"},
		SizeSelection: sizeselect.Stats{WeightedBases: 2500, WeightedFragments: 100, MeanWeightedLength: 400},
	}
	budget := SequencingBudget{ReadLayout: "pe", ReadLength: 150, LaneReadPairs: 500, Lanes: 1, UsableReadFraction: 1, Samples: 1, DesiredDepth: 10}
	target := DesignTarget{TargetGenomePct: 2.5, CoverageTolerancePct: 0.01, Objective: ObjectiveBalanced}

	cand := EvaluateSummary(summary, 100000, budget, target, DefaultScoreWeights())
	if cand.Feasible {
		t.Fatalf("candidate should not be feasible: %+v", cand)
	}
	if cand.DepthShortfallRel != 0.5 {
		t.Fatalf("depth shortfall rel = %g, want 0.5", cand.DepthShortfallRel)
	}
}

func TestSortCandidatesBalancedPrefersFeasibleThenLoss(t *testing.T) {
	candidates := []Candidate{
		{EnzymeA: "B", EnzymeB: "C", Feasible: false, DesignLoss: 0.01},
		{EnzymeA: "A", EnzymeB: "C", Feasible: true, DesignLoss: 0.2},
		{EnzymeA: "A", EnzymeB: "B", Feasible: true, DesignLoss: 0.1},
	}
	SortCandidates(candidates, ObjectiveBalanced)
	if candidates[0].EnzymeA != "A" || candidates[0].EnzymeB != "B" || candidates[0].Rank != 1 {
		t.Fatalf("unexpected top candidate: %+v", candidates[0])
	}
}

func TestCountReferenceBases(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(path, []byte(">chr1\nACGTNN\n>chr2\nNNAA\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	bases, err := CountReferenceBases(path)
	if err != nil {
		t.Fatalf("CountReferenceBases returned error: %v", err)
	}
	if bases.AllBases != 10 || bases.NonNBases != 6 {
		t.Fatalf("bases = %+v, want all=10 nonN=6", bases)
	}
}
