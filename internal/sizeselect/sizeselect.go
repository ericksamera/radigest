package sizeselect

import (
	"fmt"
	"math"
	"strings"
)

type Model string

const (
	ModelHard       Model = "hard"
	ModelNormal     Model = "normal"
	ModelTriangular Model = "triangular"
	ModelSoftWindow Model = "soft-window"
)

type Config struct {
	Model    Model   `json:"model"`
	Min      int     `json:"min_length"`
	Max      int     `json:"max_length"`
	ScoreMin int     `json:"score_min"`
	ScoreMax int     `json:"score_max"`
	Mean     float64 `json:"mean,omitempty"`
	SD       float64 `json:"sd,omitempty"`
	EdgeSD   float64 `json:"edge_sd,omitempty"`
}

type Selector struct {
	cfg Config
}

func New(cfg Config) (Selector, error) {
	cfg.Model = Model(strings.ToLower(strings.TrimSpace(string(cfg.Model))))
	if cfg.Model == "" {
		cfg.Model = ModelHard
	}
	if cfg.Min < 0 {
		return Selector{}, fmt.Errorf("-min must be >= 0 for size modeling (got %d)", cfg.Min)
	}
	if cfg.Max < cfg.Min {
		return Selector{}, fmt.Errorf("-max must be >= -min for size modeling (got min=%d max=%d)", cfg.Min, cfg.Max)
	}
	if cfg.ScoreMin < 0 {
		return Selector{}, fmt.Errorf("-score-min must be >= 0 (got %d)", cfg.ScoreMin)
	}
	if cfg.ScoreMax < cfg.ScoreMin {
		return Selector{}, fmt.Errorf("-score-max must be >= -score-min (got score-min=%d score-max=%d)", cfg.ScoreMin, cfg.ScoreMax)
	}
	if cfg.Mean == 0 {
		cfg.Mean = float64(cfg.Min+cfg.Max) / 2
	}

	switch cfg.Model {
	case ModelHard:
		// no additional parameters
	case ModelNormal:
		if !finitePositive(cfg.SD) {
			return Selector{}, fmt.Errorf("-size-sd must be > 0 for -size-model normal (got %g)", cfg.SD)
		}
		if !finite(cfg.Mean) {
			return Selector{}, fmt.Errorf("-size-mean must be finite for -size-model normal (got %g)", cfg.Mean)
		}
	case ModelTriangular:
		if !finite(cfg.Mean) || cfg.Mean <= float64(cfg.Min) || cfg.Mean >= float64(cfg.Max) {
			return Selector{}, fmt.Errorf("-size-mean must be inside (-min,-max) for -size-model triangular (got mean=%g min=%d max=%d)", cfg.Mean, cfg.Min, cfg.Max)
		}
	case ModelSoftWindow:
		if !finitePositive(cfg.EdgeSD) {
			return Selector{}, fmt.Errorf("-size-edge-sd must be > 0 for -size-model soft-window (got %g)", cfg.EdgeSD)
		}
	default:
		return Selector{}, fmt.Errorf("unknown -size-model %q; use hard, normal, triangular, or soft-window", cfg.Model)
	}

	return Selector{cfg: cfg}, nil
}

func finite(x float64) bool {
	return !math.IsNaN(x) && !math.IsInf(x, 0)
}

func finitePositive(x float64) bool {
	return finite(x) && x > 0
}

func (s Selector) Config() Config { return s.cfg }

func (s Selector) InScoreRange(length int) bool {
	return length >= s.cfg.ScoreMin && length <= s.cfg.ScoreMax
}

func (s Selector) InHardWindow(length int) bool {
	return length >= s.cfg.Min && length <= s.cfg.Max
}

func (s Selector) Weight(length int) float64 {
	l := float64(length)
	switch s.cfg.Model {
	case ModelHard:
		if s.InHardWindow(length) {
			return 1
		}
		return 0
	case ModelNormal:
		z := (l - s.cfg.Mean) / s.cfg.SD
		return math.Exp(-0.5 * z * z)
	case ModelTriangular:
		if l < float64(s.cfg.Min) || l > float64(s.cfg.Max) {
			return 0
		}
		if l == s.cfg.Mean {
			return 1
		}
		if l < s.cfg.Mean {
			return clamp01((l - float64(s.cfg.Min)) / (s.cfg.Mean - float64(s.cfg.Min)))
		}
		return clamp01((float64(s.cfg.Max) - l) / (float64(s.cfg.Max) - s.cfg.Mean))
	case ModelSoftWindow:
		left := sigmoid((l - float64(s.cfg.Min)) / s.cfg.EdgeSD)
		right := sigmoid((float64(s.cfg.Max) - l) / s.cfg.EdgeSD)
		return left * right
	default:
		return 0
	}
}

func sigmoid(x float64) float64 {
	if x >= 40 {
		return 1
	}
	if x <= -40 {
		return 0
	}
	return 1 / (1 + math.Exp(-x))
}

func clamp01(x float64) float64 {
	if x < 0 {
		return 0
	}
	if x > 1 {
		return 1
	}
	return x
}

type Stats struct {
	Model                Model   `json:"model"`
	ScoreMin             int     `json:"score_min"`
	ScoreMax             int     `json:"score_max"`
	Mean                 float64 `json:"mean,omitempty"`
	SD                   float64 `json:"sd,omitempty"`
	EdgeSD               float64 `json:"edge_sd,omitempty"`
	RawFragmentsScored   int     `json:"raw_fragments_scored"`
	RawBasesScored       int64   `json:"raw_bases_scored"`
	RawFragmentsInWindow int     `json:"raw_fragments_in_window"`
	RawBasesInWindow     int64   `json:"raw_bases_in_window"`
	WeightedFragments    float64 `json:"weighted_fragments"`
	WeightedBases        float64 `json:"weighted_bases"`
	MeanWeightedLength   float64 `json:"mean_weighted_length"`
}

func NewStats(s Selector) Stats {
	cfg := s.Config()
	st := Stats{
		Model:    cfg.Model,
		ScoreMin: cfg.ScoreMin,
		ScoreMax: cfg.ScoreMax,
	}
	switch cfg.Model {
	case ModelNormal:
		st.Mean = cfg.Mean
		st.SD = cfg.SD
	case ModelTriangular:
		st.Mean = cfg.Mean
	case ModelSoftWindow:
		st.EdgeSD = cfg.EdgeSD
	}
	return st
}

func (s *Stats) AddScored(length int, weight float64) {
	s.RawFragmentsScored++
	s.RawBasesScored += int64(length)
	s.WeightedFragments += weight
	s.WeightedBases += float64(length) * weight
	if s.WeightedFragments > 0 {
		s.MeanWeightedLength = s.WeightedBases / s.WeightedFragments
	}
}

func (s *Stats) AddHardKept(length int) {
	s.RawFragmentsInWindow++
	s.RawBasesInWindow += int64(length)
}

func (s *Stats) Add(length int, hardKept bool, weight float64) {
	s.AddScored(length, weight)
	if hardKept {
		s.AddHardKept(length)
	}
}
