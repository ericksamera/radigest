package screen

import (
	"bytes"
	"encoding/json"
	"math"
	"os"
	"path/filepath"
	"reflect"
	"testing"

	"github.com/ericksamera/radigest/internal/digest"
	"github.com/ericksamera/radigest/internal/enzyme"
	"github.com/ericksamera/radigest/internal/fasta"
	"github.com/ericksamera/radigest/internal/sizeselect"
)

func testSelector(t *testing.T) sizeselect.Selector {
	t.Helper()
	sel, err := sizeselect.New(sizeselect.Config{
		Model:    sizeselect.ModelSoftWindow,
		Min:      1,
		Max:      100,
		ScoreMin: 1,
		ScoreMax: 100,
		EdgeSD:   10,
	})
	if err != nil {
		t.Fatalf("sizeselect.New returned error: %v", err)
	}
	return sel
}

func testRecords() []fasta.Record {
	return []fasta.Record{
		{ID: "toy", Seq: []byte("AAAAGAATTCTTAAAGAATTCTTT")},
		{ID: "nocut", Seq: []byte("CCCCCCCC")},
	}
}

func testEnzymes() []enzyme.Enzyme {
	return []enzyme.Enzyme{
		{Name: "EcoRI", Recognition: "G^AATTC"},
		{Name: "MseI", Recognition: "T^TAA"},
		{Name: "ApeKI", Recognition: "G^CWGC"},
	}
}

func expectedFromPlan(t *testing.T, records []fasta.Record, ens []enzyme.Enzyme, selector sizeselect.Selector, opt digest.Options) PairSummary {
	t.Helper()
	plan, err := digest.TryNewPlanWithOptions(ens, opt)
	if err != nil {
		t.Fatalf("TryNewPlanWithOptions returned error: %v", err)
	}

	cfg := selector.Config()
	digestMin := minInt(cfg.Min, cfg.ScoreMin)
	digestMax := maxInt(cfg.Max, cfg.ScoreMax)
	sizeStats := sizeselect.NewStats(selector)
	perChromosome := make(map[string]RecordStats, len(records))
	totalFragments := 0
	totalBases := 0

	for _, rec := range records {
		local := RecordStats{}
		err := plan.DigestEach(rec.Seq, digestMin, digestMax, func(fr digest.Fragment) error {
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
			t.Fatalf("DigestEach returned error: %v", err)
		}
		perChromosome[rec.ID] = local
	}

	cfg = selector.Config()
	return PairSummary{
		Enzymes:        []string{ens[0].Name, ens[1].Name},
		MinLength:      cfg.Min,
		MaxLength:      cfg.Max,
		TotalFragments: totalFragments,
		TotalBases:     totalBases,
		PerChromosome:  perChromosome,
		SizeSelection:  sizeStats,
	}
}

func assertFloatNear(t *testing.T, name string, got, want float64) {
	t.Helper()
	if math.Abs(got-want) > 1e-9 {
		t.Fatalf("%s got %.12g want %.12g", name, got, want)
	}
}

func assertSummaryMatches(t *testing.T, got, want PairSummary) {
	t.Helper()
	if !reflect.DeepEqual(got.Enzymes, want.Enzymes) {
		t.Fatalf("enzymes got %#v want %#v", got.Enzymes, want.Enzymes)
	}
	if got.MinLength != want.MinLength || got.MaxLength != want.MaxLength {
		t.Fatalf("length window got %d-%d want %d-%d", got.MinLength, got.MaxLength, want.MinLength, want.MaxLength)
	}
	if got.TotalFragments != want.TotalFragments || got.TotalBases != want.TotalBases {
		t.Fatalf("totals got fragments=%d bases=%d want fragments=%d bases=%d", got.TotalFragments, got.TotalBases, want.TotalFragments, want.TotalBases)
	}
	if !reflect.DeepEqual(got.PerChromosome, want.PerChromosome) {
		t.Fatalf("per-chromosome stats got %#v want %#v", got.PerChromosome, want.PerChromosome)
	}
	if got.SizeSelection.RawFragmentsScored != want.SizeSelection.RawFragmentsScored ||
		got.SizeSelection.RawBasesScored != want.SizeSelection.RawBasesScored ||
		got.SizeSelection.RawFragmentsInWindow != want.SizeSelection.RawFragmentsInWindow ||
		got.SizeSelection.RawBasesInWindow != want.SizeSelection.RawBasesInWindow {
		t.Fatalf("raw size stats got %#v want %#v", got.SizeSelection, want.SizeSelection)
	}
	assertFloatNear(t, "weighted fragments", got.SizeSelection.WeightedFragments, want.SizeSelection.WeightedFragments)
	assertFloatNear(t, "weighted bases", got.SizeSelection.WeightedBases, want.SizeSelection.WeightedBases)
	assertFloatNear(t, "mean weighted length", got.SizeSelection.MeanWeightedLength, want.SizeSelection.MeanWeightedLength)
}

func TestBuildCutIndexFromRecordsMatchesBuildCutIndex(t *testing.T) {
	records := testRecords()
	ch := make(chan fasta.Record)
	go func() {
		defer close(ch)
		for _, rec := range records {
			ch <- rec
		}
	}()

	got, err := BuildCutIndexFromRecords(ch, testEnzymes(), digest.Options{})
	if err != nil {
		t.Fatalf("BuildCutIndexFromRecords returned error: %v", err)
	}
	want, err := BuildCutIndex(records, testEnzymes(), digest.Options{})
	if err != nil {
		t.Fatalf("BuildCutIndex returned error: %v", err)
	}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("streamed cut index got %#v want %#v", got, want)
	}
}

func TestBuildCutIndexFromRecordsRejectsNilChannel(t *testing.T) {
	_, err := BuildCutIndexFromRecords(nil, testEnzymes(), digest.Options{})
	if err == nil {
		t.Fatal("BuildCutIndexFromRecords accepted nil records channel")
	}
}

func TestBuildCutIndexFromFASTAMatchesBuildCutIndex(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "records.fa")
	data := []byte(`>toy
AAAAGAATTCTTAAAGAATTCTTT
>nocut
CCCCCCCC
`)
	if err := os.WriteFile(path, data, 0o644); err != nil {
		t.Fatalf("write FASTA: %v", err)
	}

	got, err := BuildCutIndexFromFASTA(path, testEnzymes(), digest.Options{})
	if err != nil {
		t.Fatalf("BuildCutIndexFromFASTA returned error: %v", err)
	}
	want, err := BuildCutIndex(testRecords(), testEnzymes(), digest.Options{})
	if err != nil {
		t.Fatalf("BuildCutIndex returned error: %v", err)
	}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("FASTA-streamed cut index got %#v want %#v", got, want)
	}
}

func TestBuildCutIndexCachesExpectedCuts(t *testing.T) {
	idx, err := BuildCutIndex(testRecords(), testEnzymes(), digest.Options{})
	if err != nil {
		t.Fatalf("BuildCutIndex returned error: %v", err)
	}
	if len(idx.Records) != 2 {
		t.Fatalf("records got %d want 2", len(idx.Records))
	}
	got := idx.Records[0].Cuts["EcoRI"]
	want := []int{5, 16}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("EcoRI cuts got %#v want %#v", got, want)
	}
	got = idx.Records[0].Cuts["MseI"]
	want = []int{11}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("MseI cuts got %#v want %#v", got, want)
	}
	if idx.CachedCutSites() != 3 {
		t.Fatalf("cached cut sites got %d want 3", idx.CachedCutSites())
	}
	if idx.CacheMemoryEstimateBytes() <= 0 {
		t.Fatalf("cache memory estimate should be positive")
	}
}

func TestBuildCutIndexRejectsDuplicateEnzymeNames(t *testing.T) {
	_, err := BuildCutIndex(testRecords(), []enzyme.Enzyme{
		{Name: "EcoRI", Recognition: "G^AATTC"},
		{Name: "EcoRI", Recognition: "G^AATTC"},
	}, digest.Options{})
	if err == nil {
		t.Fatal("BuildCutIndex accepted duplicate enzyme names")
	}
}

func TestPairNames(t *testing.T) {
	idx, err := BuildCutIndex(testRecords(), testEnzymes(), digest.Options{})
	if err != nil {
		t.Fatalf("BuildCutIndex returned error: %v", err)
	}
	want := []Pair{{A: "EcoRI", B: "MseI"}, {A: "EcoRI", B: "ApeKI"}, {A: "MseI", B: "ApeKI"}}
	if got := idx.PairNames(); !reflect.DeepEqual(got, want) {
		t.Fatalf("PairNames got %#v want %#v", got, want)
	}
}

func TestScorePairMatchesPlanDigestEach(t *testing.T) {
	records := testRecords()
	ens := testEnzymes()
	idx, err := BuildCutIndex(records, ens, digest.Options{})
	if err != nil {
		t.Fatalf("BuildCutIndex returned error: %v", err)
	}
	sel := testSelector(t)

	got, err := ScorePair(idx, "EcoRI", "MseI", sel, digest.Options{})
	if err != nil {
		t.Fatalf("ScorePair returned error: %v", err)
	}
	want := expectedFromPlan(t, records, ens[:2], sel, digest.Options{})
	assertSummaryMatches(t, got, want)
	if got.Screening.Engine != EngineCachedCutIndex {
		t.Fatalf("screening engine got %q", got.Screening.Engine)
	}
}

func TestScoreAllPairs(t *testing.T) {
	idx, err := BuildCutIndex(testRecords(), testEnzymes(), digest.Options{})
	if err != nil {
		t.Fatalf("BuildCutIndex returned error: %v", err)
	}
	summaries, err := ScoreAllPairs(idx, testSelector(t), digest.Options{})
	if err != nil {
		t.Fatalf("ScoreAllPairs returned error: %v", err)
	}
	if len(summaries) != 3 {
		t.Fatalf("ScoreAllPairs returned %d summaries, want 3", len(summaries))
	}
}

func TestScorePairMatchesPlanDigestEachWithOptions(t *testing.T) {
	records := []fasta.Record{
		{ID: "same", Seq: []byte("AAAAGAATTCAAAAGAATTCAAA")},
	}
	ens := []enzyme.Enzyme{
		{Name: "EcoRI", Recognition: "G^AATTC"},
		{Name: "NcoI", Recognition: "C^CATGG"},
	}
	opt := digest.Options{AllowSame: true, IncludeEnds: true}
	idx, err := BuildCutIndex(records, ens, opt)
	if err != nil {
		t.Fatalf("BuildCutIndex returned error: %v", err)
	}
	sel := testSelector(t)

	got, err := ScorePair(idx, "EcoRI", "NcoI", sel, opt)
	if err != nil {
		t.Fatalf("ScorePair returned error: %v", err)
	}
	want := expectedFromPlan(t, records, ens, sel, opt)
	assertSummaryMatches(t, got, want)
}

func TestScorePairMatchesPlanDigestEachWithCoincidentCuts(t *testing.T) {
	records := []fasta.Record{{ID: "coincident", Seq: []byte("AAAGATCAAAGATC")}}
	ens := []enzyme.Enzyme{
		{Name: "DpnII", Recognition: "^GATC"},
		{Name: "MboI", Recognition: "^GATC"},
	}
	idx, err := BuildCutIndex(records, ens, digest.Options{})
	if err != nil {
		t.Fatalf("BuildCutIndex returned error: %v", err)
	}
	sel, err := sizeselect.New(sizeselect.Config{Model: sizeselect.ModelHard, Min: 0, Max: 100, ScoreMin: 0, ScoreMax: 100})
	if err != nil {
		t.Fatalf("sizeselect.New returned error: %v", err)
	}
	got, err := ScorePair(idx, "DpnII", "MboI", sel, digest.Options{})
	if err != nil {
		t.Fatalf("ScorePair returned error: %v", err)
	}
	want := expectedFromPlan(t, records, ens, sel, digest.Options{})
	assertSummaryMatches(t, got, want)
}

func TestScorePairReportsMissingEnzyme(t *testing.T) {
	idx, err := BuildCutIndex(testRecords(), testEnzymes(), digest.Options{})
	if err != nil {
		t.Fatalf("BuildCutIndex returned error: %v", err)
	}
	_, err = ScorePair(idx, "EcoRI", "Missing", testSelector(t), digest.Options{})
	if err == nil {
		t.Fatal("ScorePair accepted missing enzyme")
	}
}

func TestWritePairSummaryJSON(t *testing.T) {
	idx, err := BuildCutIndex(testRecords(), testEnzymes(), digest.Options{})
	if err != nil {
		t.Fatalf("BuildCutIndex returned error: %v", err)
	}
	summary, err := ScorePair(idx, "EcoRI", "MseI", testSelector(t), digest.Options{})
	if err != nil {
		t.Fatalf("ScorePair returned error: %v", err)
	}

	var buf bytes.Buffer
	if err := WritePairSummaryJSON(&buf, summary); err != nil {
		t.Fatalf("WritePairSummaryJSON returned error: %v", err)
	}

	var decoded struct {
		Enzymes       []string         `json:"enzymes"`
		SizeSelection sizeselect.Stats `json:"size_selection"`
		Screening     ScreeningStats   `json:"screening"`
	}
	if err := json.Unmarshal(buf.Bytes(), &decoded); err != nil {
		t.Fatalf("summary JSON did not unmarshal: %v", err)
	}
	if !reflect.DeepEqual(decoded.Enzymes, []string{"EcoRI", "MseI"}) {
		t.Fatalf("decoded enzymes got %#v", decoded.Enzymes)
	}
	if decoded.SizeSelection.RawFragmentsInWindow != summary.SizeSelection.RawFragmentsInWindow {
		t.Fatalf("decoded size selection mismatch")
	}
	if decoded.Screening.Engine != EngineCachedCutIndex {
		t.Fatalf("decoded engine got %q", decoded.Screening.Engine)
	}
}
