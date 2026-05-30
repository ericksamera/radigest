// Package screen provides cached cut-index helpers for enzyme-pair screening.
//
// The package is intentionally internal and does not define a command-line
// interface. It supports cached screening in three explicit steps:
//
//  1. Scan each candidate enzyme once per FASTA record and store cut positions.
//  2. Reuse those cut streams for many enzyme-pair scoring operations.
//  3. Produce pair summaries with the core fields consumed by
//     scripts/radigest-rank-pairs.
package screen

import (
	"encoding/json"
	"fmt"
	"io"
	"strconv"

	"github.com/ericksamera/radigest/internal/digest"
	"github.com/ericksamera/radigest/internal/enzyme"
	"github.com/ericksamera/radigest/internal/fasta"
	"github.com/ericksamera/radigest/internal/sizeselect"
)

const EngineCachedCutIndex = "cached-cut-index"

// RecordCuts stores cached cut coordinates for one FASTA record.
type RecordCuts struct {
	ID     string
	Length int
	Cuts   map[string][]int
}

// CutIndex stores per-record, per-enzyme sorted cut coordinates.
type CutIndex struct {
	Records     []RecordCuts
	EnzymeNames []string
}

// RecordStats summarizes hard-window fragments for one record.
type RecordStats struct {
	Fragments int `json:"fragments"`
	Bases     int `json:"bases"`
}

// ScreeningStats records cached-screening provenance and cache size information.
type ScreeningStats struct {
	Engine                   string `json:"engine"`
	CandidateEnzymes         int    `json:"candidate_enzymes"`
	Records                  int    `json:"records"`
	CachedCutSites           int    `json:"cached_cut_sites"`
	CacheMemoryEstimateBytes int64  `json:"cache_memory_estimate_bytes"`
}

// PairSummary is the JSON-compatible summary emitted by cached pair screening.
// It intentionally preserves fields consumed by scripts/radigest-rank-pairs:
// enzymes and size_selection.raw/weighted fields.
type PairSummary struct {
	SchemaVersion  int                    `json:"schema_version"`
	Enzymes        []string               `json:"enzymes"`
	MinLength      int                    `json:"min_length"`
	MaxLength      int                    `json:"max_length"`
	TotalFragments int                    `json:"total_fragments"`
	TotalBases     int                    `json:"total_bases"`
	PerChromosome  map[string]RecordStats `json:"per_chromosome"`
	SizeSelection  sizeselect.Stats       `json:"size_selection"`
	Screening      ScreeningStats         `json:"screening"`
}

// Pair identifies one unique enzyme pair.
type Pair struct {
	A string
	B string
}

// BuildCutIndex scans every input record once per enzyme and stores sorted cut
// coordinates. Enzyme names must be unique because they are used as map keys.
//
// BuildCutIndex accepts a materialized record slice for tests and callers that
// already hold FASTA records in memory. Cached screening CLIs should prefer
// BuildCutIndexFromFASTA or BuildCutIndexFromRecords so sequence bytes can be
// discarded after each record is scanned.
func BuildCutIndex(records []fasta.Record, enzymes []enzyme.Enzyme, opt digest.Options) (CutIndex, error) {
	names, plans, err := compileCutPlans(enzymes, opt)
	if err != nil {
		return CutIndex{}, err
	}

	idx := CutIndex{
		Records:     make([]RecordCuts, 0, len(records)),
		EnzymeNames: names,
	}

	for _, rec := range records {
		rc, err := scanRecordCuts(rec, names, plans)
		if err != nil {
			return CutIndex{}, err
		}
		idx.Records = append(idx.Records, rc)
	}

	return idx, nil
}

// BuildCutIndexFromRecords scans records from a channel and stores only record
// IDs, record lengths, and per-enzyme cut coordinates. Sequence bytes are not
// retained after each record has been scanned.
func BuildCutIndexFromRecords(records <-chan fasta.Record, enzymes []enzyme.Enzyme, opt digest.Options) (CutIndex, error) {
	if records == nil {
		return CutIndex{}, fmt.Errorf("screen cut index: records channel is nil")
	}

	names, plans, err := compileCutPlans(enzymes, opt)
	if err != nil {
		return CutIndex{}, err
	}

	idx := CutIndex{
		Records:     make([]RecordCuts, 0),
		EnzymeNames: names,
	}

	for rec := range records {
		rc, err := scanRecordCuts(rec, names, plans)
		if err != nil {
			return CutIndex{}, err
		}
		idx.Records = append(idx.Records, rc)
	}

	return idx, nil
}

// BuildCutIndexFromFASTA streams a FASTA file and builds a cut index without
// retaining full reference sequences. This is the preferred entry point for
// cached pair-screening commands on large genomes.
func BuildCutIndexFromFASTA(path string, enzymes []enzyme.Enzyme, opt digest.Options) (CutIndex, error) {
	ch := make(chan fasta.Record)
	errCh := make(chan error, 1)
	go func() {
		errCh <- fasta.Stream(path, ch)
	}()

	idx, buildErr := BuildCutIndexFromRecords(ch, enzymes, opt)
	if buildErr != nil {
		for range ch {
			// Drain the FASTA stream so fasta.Stream can return its error and the
			// goroutine cannot block while sending a later record.
		}
	}
	streamErr := <-errCh
	if buildErr != nil {
		return CutIndex{}, buildErr
	}
	if streamErr != nil {
		return CutIndex{}, streamErr
	}
	return idx, nil
}

func compileCutPlans(enzymes []enzyme.Enzyme, opt digest.Options) ([]string, []digest.Plan, error) {
	if len(enzymes) == 0 {
		return nil, nil, fmt.Errorf("screen cut index: no enzymes provided")
	}

	names := make([]string, 0, len(enzymes))
	plans := make([]digest.Plan, 0, len(enzymes))
	seen := make(map[string]struct{}, len(enzymes))

	for _, enz := range enzymes {
		if enz.Name == "" {
			return nil, nil, fmt.Errorf("screen cut index: enzyme with empty name")
		}
		if _, ok := seen[enz.Name]; ok {
			return nil, nil, fmt.Errorf("screen cut index: duplicate enzyme name %q", enz.Name)
		}
		seen[enz.Name] = struct{}{}
		names = append(names, enz.Name)

		plan, err := digest.TryNewPlanWithOptions([]enzyme.Enzyme{enz}, digest.Options{StrictCuts: opt.StrictCuts})
		if err != nil {
			return nil, nil, err
		}
		plans = append(plans, plan)
	}

	return names, plans, nil
}

func scanRecordCuts(rec fasta.Record, names []string, plans []digest.Plan) (RecordCuts, error) {
	if rec.ID == "" {
		return RecordCuts{}, fmt.Errorf("screen cut index: record with empty ID")
	}
	rc := RecordCuts{
		ID:     rec.ID,
		Length: len(rec.Seq),
		Cuts:   make(map[string][]int, len(names)),
	}
	for i, plan := range plans {
		rc.Cuts[names[i]] = plan.Cuts(rec.Seq)
	}
	return rc, nil
}

// ContainsEnzyme reports whether the cut index includes name.
func (idx CutIndex) ContainsEnzyme(name string) bool {
	for _, enzymeName := range idx.EnzymeNames {
		if enzymeName == name {
			return true
		}
	}
	return false
}

// PairNames returns all unique enzyme pairs in cut-index enzyme order.
func (idx CutIndex) PairNames() []Pair {
	pairs := make([]Pair, 0, len(idx.EnzymeNames)*(len(idx.EnzymeNames)-1)/2)
	for i := 0; i < len(idx.EnzymeNames); i++ {
		for j := i + 1; j < len(idx.EnzymeNames); j++ {
			pairs = append(pairs, Pair{A: idx.EnzymeNames[i], B: idx.EnzymeNames[j]})
		}
	}
	return pairs
}

// CachedCutSites returns the total number of cached cut coordinates.
func (idx CutIndex) CachedCutSites() int {
	total := 0
	for _, rec := range idx.Records {
		for _, cuts := range rec.Cuts {
			total += len(cuts)
		}
	}
	return total
}

// CacheMemoryEstimateBytes returns an approximate in-memory size for cached cut
// coordinates only. It intentionally excludes map, slice, and string overhead.
func (idx CutIndex) CacheMemoryEstimateBytes() int64 {
	return int64(idx.CachedCutSites()) * int64(strconv.IntSize/8)
}

// ScorePair scores one enzyme pair from cached cut-coordinate streams.
func ScorePair(idx CutIndex, enzymeA, enzymeB string, selector sizeselect.Selector, opt digest.Options) (PairSummary, error) {
	if enzymeA == "" || enzymeB == "" {
		return PairSummary{}, fmt.Errorf("screen score pair: enzyme names must be non-empty")
	}
	if enzymeA == enzymeB {
		return PairSummary{}, fmt.Errorf("screen score pair: self-pair %q is not supported", enzymeA)
	}
	if !idx.ContainsEnzyme(enzymeA) {
		return PairSummary{}, fmt.Errorf("screen score pair: enzyme %q not found in cut index", enzymeA)
	}
	if !idx.ContainsEnzyme(enzymeB) {
		return PairSummary{}, fmt.Errorf("screen score pair: enzyme %q not found in cut index", enzymeB)
	}

	cfg := selector.Config()
	digestMin := minInt(cfg.Min, cfg.ScoreMin)
	digestMax := maxInt(cfg.Max, cfg.ScoreMax)

	sizeStats := sizeselect.NewStats(selector)
	perChromosome := make(map[string]RecordStats, len(idx.Records))
	totalFragments := 0
	totalBases := 0

	for _, rec := range idx.Records {
		cutsA := rec.Cuts[enzymeA]
		cutsB := rec.Cuts[enzymeB]
		local := RecordStats{}

		err := digest.DigestCutsEach(cutsA, cutsB, rec.Length, digestMin, digestMax, opt, func(fr digest.Fragment) error {
			length := fr.End - fr.Start
			hardKept := selector.InHardWindow(length)
			if hardKept {
				sizeStats.AddHardKept(length)
				local.Fragments++
				local.Bases += length
				totalFragments++
				totalBases += length
			}
			if selector.InScoreRange(length) {
				sizeStats.AddScored(length, selector.Weight(length))
			}
			return nil
		})
		if err != nil {
			return PairSummary{}, err
		}
		perChromosome[rec.ID] = local
	}

	return PairSummary{
		SchemaVersion:  1,
		Enzymes:        []string{enzymeA, enzymeB},
		MinLength:      cfg.Min,
		MaxLength:      cfg.Max,
		TotalFragments: totalFragments,
		TotalBases:     totalBases,
		PerChromosome:  perChromosome,
		SizeSelection:  sizeStats,
		Screening: ScreeningStats{
			Engine:                   EngineCachedCutIndex,
			CandidateEnzymes:         len(idx.EnzymeNames),
			Records:                  len(idx.Records),
			CachedCutSites:           idx.CachedCutSites(),
			CacheMemoryEstimateBytes: idx.CacheMemoryEstimateBytes(),
		},
	}, nil
}

// ScoreAllPairs scores all unique enzyme pairs in cut-index order.
func ScoreAllPairs(idx CutIndex, selector sizeselect.Selector, opt digest.Options) ([]PairSummary, error) {
	pairs := idx.PairNames()
	summaries := make([]PairSummary, 0, len(pairs))
	for _, pair := range pairs {
		summary, err := ScorePair(idx, pair.A, pair.B, selector, opt)
		if err != nil {
			return nil, err
		}
		summaries = append(summaries, summary)
	}
	return summaries, nil
}

// WritePairSummaryJSON writes one pair summary as indented JSON.
func WritePairSummaryJSON(w io.Writer, summary PairSummary) error {
	if w == nil {
		return fmt.Errorf("screen pair summary writer is nil")
	}
	enc := json.NewEncoder(w)
	enc.SetIndent("", "  ")
	return enc.Encode(summary)
}

func minInt(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func maxInt(a, b int) int {
	if a > b {
		return a
	}
	return b
}
