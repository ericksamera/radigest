// Package design provides inverse enzyme-pair design calculations layered on
// top of cached radigest pair-screen summaries.
package design

import (
	"fmt"
	"math"
	"sort"
	"strings"

	"github.com/ericksamera/radigest/internal/fasta"
	"github.com/ericksamera/radigest/internal/screen"
)

const SchemaVersion = 1

type Objective string

const (
	ObjectiveBalanced               Objective = "balanced"
	ObjectiveClosestCoverage        Objective = "closest-coverage"
	ObjectiveDepthFirst             Objective = "depth-first"
	ObjectiveFeasibleLowestCoverage Objective = "feasible-lowest-coverage"
	ObjectiveMaxDepth               Objective = "max-depth"
)

// GenomeBases records reference-size denominators used for genome-percentage
// calculations. NonNBases follows the existing helper-script convention: every
// non-N FASTA character contributes to the denominator.
type GenomeBases struct {
	AllBases  int64 `json:"all_bases"`
	NonNBases int64 `json:"non_n_bases"`
}

type SequencingBudget struct {
	ReadLayout           string  `json:"read_layout"`
	ReadLength           int     `json:"read_length"`
	LaneReadPairs        float64 `json:"lane_read_pairs"`
	Lanes                int     `json:"lanes"`
	UsableReadFraction   float64 `json:"usable_read_fraction"`
	Samples              int     `json:"samples"`
	TargetMeanLocusDepth float64 `json:"target_mean_locus_depth"`
}

func (b SequencingBudget) EffectiveReadPairsPerLane() float64 {
	return b.LaneReadPairs * b.UsableReadFraction
}

func (b SequencingBudget) EffectiveReadPairsTotal() float64 {
	return b.EffectiveReadPairsPerLane() * float64(b.Lanes)
}

func (b SequencingBudget) ReadPairsPerSample() float64 {
	if b.Samples <= 0 {
		return 0
	}
	return b.EffectiveReadPairsTotal() / float64(b.Samples)
}

type DesignTarget struct {
	TargetGenomePct      float64   `json:"target_genome_pct"`
	CoverageTolerancePct float64   `json:"coverage_tolerance_pct"`
	Objective            Objective `json:"objective"`
}

type ScoreWeights struct {
	Coverage     float64 `json:"coverage"`
	Depth        float64 `json:"depth"`
	Overcoverage float64 `json:"overcoverage"`
	Insert       float64 `json:"insert"`
}

func DefaultScoreWeights() ScoreWeights {
	return ScoreWeights{
		Coverage:     1.0,
		Depth:        2.0,
		Overcoverage: 0.5,
		Insert:       0.25,
	}
}

type Candidate struct {
	Rank           int      `json:"rank"`
	EnzymeA        string   `json:"enzyme_a"`
	EnzymeB        string   `json:"enzyme_b"`
	Enzymes        []string `json:"enzymes"`
	Feasible       bool     `json:"feasible"`
	DecisionReason string   `json:"decision_reason"`

	FitScore float64 `json:"fit_score"`
	FitLoss  float64 `json:"fit_loss"`

	TargetGenomePct              float64 `json:"target_genome_pct"`
	PredictedWeightedGenomePct   float64 `json:"predicted_weighted_genome_pct"`
	CoverageErrorPctPoints       float64 `json:"coverage_error_pct_points"`
	CoverageErrorRel             float64 `json:"coverage_error_rel"`
	OvercoverageRel              float64 `json:"overcoverage_rel"`
	UndercoverageRel             float64 `json:"undercoverage_rel"`
	TargetMeanLocusDepth         float64 `json:"target_mean_locus_depth"`
	PredictedMeanLocusDepth      float64 `json:"predicted_mean_locus_depth"`
	DepthMargin                  float64 `json:"depth_margin"`
	DepthShortfallRel            float64 `json:"depth_shortfall_rel"`
	ReadPairsPerSample           float64 `json:"read_pairs_per_sample"`
	RequiredPairsPerSampleTarget float64 `json:"required_pairs_per_sample_full_target"`

	WeightedBases        float64 `json:"weighted_bases"`
	WeightedFragments    float64 `json:"weighted_fragments"`
	MeanWeightedLength   float64 `json:"mean_weighted_length"`
	RawBasesInWindow     int64   `json:"raw_bases_in_window"`
	RawFragmentsInWindow int     `json:"raw_fragments_in_window"`

	BudgetSupportedGenomePct     float64 `json:"budget_supported_genome_pct"`
	BudgetSupportedWeightedBases float64 `json:"budget_supported_weighted_bases"`
	MaxSamplesPerLaneFullTarget  int     `json:"max_samples_per_lane_full_target"`
	MaxSamplesTotalFullTarget    int     `json:"max_samples_total_full_target"`
	LanesRequiredFullTarget      int     `json:"lanes_required_full_target"`
	AdapterThresholdBP           int     `json:"adapter_threshold_bp"`
	OverlapThresholdBP           int     `json:"overlap_threshold_bp,omitempty"`
	MeanInsertCategory           string  `json:"mean_insert_category"`
	InsertPenalty                float64 `json:"insert_penalty"`
	Records                      int     `json:"records"`
	CachedCutSites               int     `json:"cached_cut_sites"`
	CacheMemoryEstimateBytes     int64   `json:"cache_memory_estimate_bytes"`
}

func CountReferenceBases(path string) (GenomeBases, error) {
	records := make(chan fasta.Record)
	errCh := make(chan error, 1)
	go func() {
		errCh <- fasta.Stream(path, records)
	}()

	var bases GenomeBases
	for rec := range records {
		bases.AllBases += int64(len(rec.Seq))
		for _, b := range rec.Seq {
			if b != 'N' {
				bases.NonNBases++
			}
		}
	}
	if err := <-errCh; err != nil {
		return GenomeBases{}, err
	}
	if bases.AllBases == 0 {
		return GenomeBases{}, fmt.Errorf("no FASTA bases found in %q", path)
	}
	return bases, nil
}

func ValidateObjective(value string) (Objective, error) {
	obj := Objective(strings.ToLower(strings.TrimSpace(value)))
	switch obj {
	case ObjectiveBalanced, ObjectiveClosestCoverage, ObjectiveDepthFirst, ObjectiveFeasibleLowestCoverage, ObjectiveMaxDepth:
		return obj, nil
	default:
		return "", fmt.Errorf("invalid objective %q; use balanced, closest-coverage, depth-first, feasible-lowest-coverage, or max-depth", value)
	}
}

func EvaluateSummary(summary screen.PairSummary, genomeBases int64, budget SequencingBudget, target DesignTarget, weights ScoreWeights) Candidate {
	weightedBases := summary.SizeSelection.WeightedBases
	weightedFragments := summary.SizeSelection.WeightedFragments
	weightedGenomePct := 0.0
	if genomeBases > 0 {
		weightedGenomePct = 100.0 * weightedBases / float64(genomeBases)
	}

	readPairsPerSample := budget.ReadPairsPerSample()
	expectedDepth := safeDiv(readPairsPerSample, weightedFragments)
	requiredPairsPerSample := budget.TargetMeanLocusDepth * weightedFragments

	coverageDelta := weightedGenomePct - target.TargetGenomePct
	coverageErrorPctPoints := math.Abs(coverageDelta)
	coverageErrorRel := safeDiv(coverageErrorPctPoints, target.TargetGenomePct)
	overcoverageRel := safeDiv(math.Max(0, coverageDelta), target.TargetGenomePct)
	undercoverageRel := safeDiv(math.Max(0, -coverageDelta), target.TargetGenomePct)
	depthMargin := expectedDepth - budget.TargetMeanLocusDepth
	depthShortfallRel := safeDiv(math.Max(0, -depthMargin), budget.TargetMeanLocusDepth)

	adapterThreshold, overlapThreshold, meanInsertCategory := MeanInsertCategory(budget.ReadLayout, budget.ReadLength, summary.SizeSelection.MeanWeightedLength)
	insertPenalty := InsertPenalty(meanInsertCategory)

	budgetSupportedGenomePct := 0.0
	budgetSupportedWeightedBases := 0.0
	if budget.TargetMeanLocusDepth > 0 && expectedDepth > 0 {
		supportedFraction := math.Min(1.0, expectedDepth/budget.TargetMeanLocusDepth)
		budgetSupportedGenomePct = weightedGenomePct * supportedFraction
		budgetSupportedWeightedBases = weightedBases * supportedFraction
	}

	maxSamplesPerLane := 0
	maxSamplesTotal := 0
	if requiredPairsPerSample > 0 {
		maxSamplesPerLane = floorNonnegative(budget.EffectiveReadPairsPerLane() / requiredPairsPerSample)
		maxSamplesTotal = floorNonnegative(budget.EffectiveReadPairsTotal() / requiredPairsPerSample)
	}
	lanesRequired := 0
	if budget.Samples > 0 && requiredPairsPerSample > 0 && budget.EffectiveReadPairsPerLane() > 0 {
		lanesRequired = ceilNonnegative(float64(budget.Samples) * requiredPairsPerSample / budget.EffectiveReadPairsPerLane())
	}

	coverageCloseEnough := coverageErrorPctPoints <= target.CoverageTolerancePct
	depthSufficient := expectedDepth >= budget.TargetMeanLocusDepth
	feasible := weightedFragments > 0 && coverageCloseEnough && depthSufficient

	loss := weights.Coverage*coverageErrorRel +
		weights.Depth*depthShortfallRel +
		weights.Overcoverage*overcoverageRel +
		weights.Insert*insertPenalty

	score := 0.0
	if loss >= 0 && !math.IsNaN(loss) && !math.IsInf(loss, 0) {
		score = 1.0 / (1.0 + loss)
	}

	candidate := Candidate{
		Enzymes:                      append([]string(nil), summary.Enzymes...),
		Feasible:                     feasible,
		FitScore:                     score,
		FitLoss:                      loss,
		TargetGenomePct:              target.TargetGenomePct,
		PredictedWeightedGenomePct:   weightedGenomePct,
		CoverageErrorPctPoints:       coverageErrorPctPoints,
		CoverageErrorRel:             coverageErrorRel,
		OvercoverageRel:              overcoverageRel,
		UndercoverageRel:             undercoverageRel,
		TargetMeanLocusDepth:         budget.TargetMeanLocusDepth,
		PredictedMeanLocusDepth:      expectedDepth,
		DepthMargin:                  depthMargin,
		DepthShortfallRel:            depthShortfallRel,
		ReadPairsPerSample:           readPairsPerSample,
		RequiredPairsPerSampleTarget: requiredPairsPerSample,
		WeightedBases:                weightedBases,
		WeightedFragments:            weightedFragments,
		MeanWeightedLength:           summary.SizeSelection.MeanWeightedLength,
		RawBasesInWindow:             summary.SizeSelection.RawBasesInWindow,
		RawFragmentsInWindow:         summary.SizeSelection.RawFragmentsInWindow,
		BudgetSupportedGenomePct:     budgetSupportedGenomePct,
		BudgetSupportedWeightedBases: budgetSupportedWeightedBases,
		MaxSamplesPerLaneFullTarget:  maxSamplesPerLane,
		MaxSamplesTotalFullTarget:    maxSamplesTotal,
		LanesRequiredFullTarget:      lanesRequired,
		AdapterThresholdBP:           adapterThreshold,
		OverlapThresholdBP:           overlapThreshold,
		MeanInsertCategory:           meanInsertCategory,
		InsertPenalty:                insertPenalty,
		Records:                      summary.Screening.Records,
		CachedCutSites:               summary.Screening.CachedCutSites,
		CacheMemoryEstimateBytes:     summary.Screening.CacheMemoryEstimateBytes,
	}
	if len(summary.Enzymes) > 0 {
		candidate.EnzymeA = summary.Enzymes[0]
	}
	if len(summary.Enzymes) > 1 {
		candidate.EnzymeB = summary.Enzymes[1]
	}
	candidate.DecisionReason = DecisionReason(candidate, target)
	return candidate
}

func MeanInsertCategory(readLayout string, readLength int, meanWeightedLength float64) (int, int, string) {
	adapterThreshold := readLength
	if strings.ToLower(strings.TrimSpace(readLayout)) == "pe" {
		overlapThreshold := 2 * readLength
		switch {
		case meanWeightedLength <= 0:
			return adapterThreshold, overlapThreshold, "unknown"
		case meanWeightedLength < float64(adapterThreshold):
			return adapterThreshold, overlapThreshold, "mean_lt_read_length_adapter_risk"
		case meanWeightedLength < float64(overlapThreshold):
			return adapterThreshold, overlapThreshold, "mean_lt_2_read_lengths_overlap_risk"
		default:
			return adapterThreshold, overlapThreshold, "mean_ge_2_read_lengths"
		}
	}

	switch {
	case meanWeightedLength <= 0:
		return adapterThreshold, 0, "unknown"
	case meanWeightedLength < float64(adapterThreshold):
		return adapterThreshold, 0, "mean_lt_read_length_adapter_risk"
	default:
		return adapterThreshold, 0, "mean_ge_read_length"
	}
}

func InsertPenalty(category string) float64 {
	switch category {
	case "mean_ge_2_read_lengths", "mean_ge_read_length":
		return 0
	case "mean_lt_2_read_lengths_overlap_risk":
		return 0.25
	case "mean_lt_read_length_adapter_risk":
		return 1.0
	default:
		return 0.5
	}
}

func DecisionReason(c Candidate, target DesignTarget) string {
	if c.WeightedFragments <= 0 {
		return "no recovered weighted fragments"
	}
	parts := make([]string, 0, 4)
	if c.CoverageErrorPctPoints <= target.CoverageTolerancePct {
		parts = append(parts, "matches target coverage")
	} else if c.PredictedWeightedGenomePct < c.TargetGenomePct {
		parts = append(parts, fmt.Sprintf("under target coverage by %.6f pct-points", c.TargetGenomePct-c.PredictedWeightedGenomePct))
	} else {
		parts = append(parts, fmt.Sprintf("over target coverage by %.6f pct-points", c.PredictedWeightedGenomePct-c.TargetGenomePct))
	}
	if c.PredictedMeanLocusDepth >= c.TargetMeanLocusDepth {
		parts = append(parts, "meets target mean locus depth")
	} else {
		parts = append(parts, fmt.Sprintf("depth shortfall %.6f", c.TargetMeanLocusDepth-c.PredictedMeanLocusDepth))
	}
	if c.MeanInsertCategory == "mean_lt_read_length_adapter_risk" || c.MeanInsertCategory == "mean_lt_2_read_lengths_overlap_risk" {
		parts = append(parts, c.MeanInsertCategory)
	}
	return strings.Join(parts, "; ")
}

func SortCandidates(candidates []Candidate, objective Objective) {
	sort.SliceStable(candidates, func(i, j int) bool {
		a := candidates[i]
		b := candidates[j]
		switch objective {
		case ObjectiveClosestCoverage:
			if a.CoverageErrorRel != b.CoverageErrorRel {
				return a.CoverageErrorRel < b.CoverageErrorRel
			}
			if a.DepthShortfallRel != b.DepthShortfallRel {
				return a.DepthShortfallRel < b.DepthShortfallRel
			}
		case ObjectiveDepthFirst:
			if a.DepthShortfallRel != b.DepthShortfallRel {
				return a.DepthShortfallRel < b.DepthShortfallRel
			}
			if a.CoverageErrorRel != b.CoverageErrorRel {
				return a.CoverageErrorRel < b.CoverageErrorRel
			}
		case ObjectiveFeasibleLowestCoverage:
			if a.Feasible != b.Feasible {
				return a.Feasible
			}
			if a.PredictedWeightedGenomePct != b.PredictedWeightedGenomePct {
				return a.PredictedWeightedGenomePct < b.PredictedWeightedGenomePct
			}
		case ObjectiveMaxDepth:
			if a.PredictedMeanLocusDepth != b.PredictedMeanLocusDepth {
				return a.PredictedMeanLocusDepth > b.PredictedMeanLocusDepth
			}
			if a.CoverageErrorRel != b.CoverageErrorRel {
				return a.CoverageErrorRel < b.CoverageErrorRel
			}
		default:
			if a.Feasible != b.Feasible {
				return a.Feasible
			}
			if a.FitLoss != b.FitLoss {
				return a.FitLoss < b.FitLoss
			}
			if a.CoverageErrorRel != b.CoverageErrorRel {
				return a.CoverageErrorRel < b.CoverageErrorRel
			}
		}
		if a.FitLoss != b.FitLoss {
			return a.FitLoss < b.FitLoss
		}
		if a.EnzymeA != b.EnzymeA {
			return a.EnzymeA < b.EnzymeA
		}
		return a.EnzymeB < b.EnzymeB
	})
	for i := range candidates {
		candidates[i].Rank = i + 1
	}
}

func safeDiv(numerator, denominator float64) float64 {
	if denominator <= 0 || math.IsNaN(denominator) || math.IsInf(denominator, 0) {
		return 0
	}
	return numerator / denominator
}

func floorNonnegative(value float64) int {
	if value < 0 || math.IsNaN(value) || math.IsInf(value, 0) {
		return 0
	}
	return int(math.Floor(value))
}

func ceilNonnegative(value float64) int {
	if value < 0 || math.IsNaN(value) || math.IsInf(value, 0) {
		return 0
	}
	return int(math.Ceil(value))
}
