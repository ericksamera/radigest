package main

import (
	"encoding/csv"
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"io"
	"math"
	"os"
	"path/filepath"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"sync"

	"github.com/ericksamera/radigest/internal/design"
	"github.com/ericksamera/radigest/internal/digest"
	"github.com/ericksamera/radigest/internal/enzyme"
	"github.com/ericksamera/radigest/internal/screen"
	"github.com/ericksamera/radigest/internal/sizeselect"
)

var version = "dev"

type usageError struct{ err error }

func (e usageError) Error() string { return e.err.Error() }
func (e usageError) Unwrap() error { return e.err }

type cliConfig struct {
	fastaPath            string
	enzFlag              string
	outDir               string
	tsvPath              string
	jsonPath             string
	denominator          string
	genomeBases          int64
	minLen               int
	maxLen               int
	scoreMin             int
	scoreMax             int
	sizeModel            string
	sizeMean             float64
	sizeSD               float64
	sizeEdgeSD           float64
	allowSame            bool
	includeEnds          bool
	strictCuts           bool
	readLayout           string
	readLength           int
	laneReadPairs        float64
	flowcellReadPairs    float64
	lanes                int
	usableReadFraction   float64
	samples              int
	desiredDepth         float64
	targetGenomePct      float64
	coverageTolerancePct float64
	objective            string
	weightCoverage       float64
	weightDepth          float64
	weightOvercoverage   float64
	weightInsert         float64
	jobs                 int
	threads              int
	buildWorkers         int
	maxPairs             int
	top                  int
	force                bool
	showVersion          bool
	help                 bool
}

type digestParameters struct {
	MinLength   int     `json:"min_length"`
	MaxLength   int     `json:"max_length"`
	ScoreMin    int     `json:"score_min"`
	ScoreMax    int     `json:"score_max"`
	SizeModel   string  `json:"size_model"`
	SizeMean    float64 `json:"size_mean,omitempty"`
	SizeSD      float64 `json:"size_sd,omitempty"`
	SizeEdgeSD  float64 `json:"size_edge_sd,omitempty"`
	AllowSame   bool    `json:"allow_same"`
	IncludeEnds bool    `json:"include_ends"`
	StrictCuts  bool    `json:"strict_cuts"`
}

type inputSummary struct {
	FASTA       string             `json:"fasta"`
	Denominator string             `json:"denominator"`
	GenomeBases int64              `json:"genome_bases"`
	Reference   design.GenomeBases `json:"reference_bases"`
}

type runSummary struct {
	CandidateEnzymes int      `json:"candidate_enzymes"`
	CandidatePairs   int      `json:"candidate_pairs"`
	ReportedPairs    int      `json:"reported_pairs"`
	FeasiblePairs    int      `json:"feasible_pairs"`
	BestPair         []string `json:"best_pair,omitempty"`
}

type outputSummary struct {
	TSV  string `json:"tsv"`
	JSON string `json:"json"`
}

type designReport struct {
	SchemaVersion   int                     `json:"schema_version"`
	RadigestVersion string                  `json:"radigest_version"`
	Command         []string                `json:"command"`
	Input           inputSummary            `json:"input"`
	Digest          digestParameters        `json:"digest_parameters"`
	Sequencing      design.SequencingBudget `json:"sequencing_budget"`
	Target          design.DesignTarget     `json:"design_target"`
	Weights         design.ScoreWeights     `json:"score_weights"`
	Screening       screen.ScreeningStats   `json:"screening"`
	Outputs         outputSummary           `json:"outputs"`
	Warnings        []string                `json:"warnings"`
	Summary         runSummary              `json:"summary"`
	Results         []design.Candidate      `json:"results"`
}

func main() {
	if err := run(os.Args[1:], os.Stdout, os.Stderr); err != nil {
		code := 1
		var usage usageError
		if errors.As(err, &usage) {
			code = 2
		}
		_, _ = fmt.Fprintln(os.Stderr, "error:", err)
		os.Exit(code)
	}
}

func run(args []string, stdout, stderr io.Writer) error {
	if stdout == nil {
		stdout = io.Discard
	}
	if stderr == nil {
		stderr = io.Discard
	}

	cfg, err := parseArgs(args, stderr)
	if err != nil {
		return err
	}
	if cfg.help {
		return nil
	}
	if cfg.showVersion {
		_, err := fmt.Fprintf(stdout, "radigest-design %s\n", version)
		return err
	}

	objective, err := design.ValidateObjective(cfg.objective)
	if err != nil {
		return usageError{err: err}
	}
	selector, err := sizeselect.New(sizeselect.Config{
		Model:    sizeselect.Model(cfg.sizeModel),
		Min:      cfg.minLen,
		Max:      cfg.maxLen,
		ScoreMin: cfg.scoreMin,
		ScoreMax: cfg.scoreMax,
		Mean:     cfg.sizeMean,
		SD:       cfg.sizeSD,
		EdgeSD:   cfg.sizeEdgeSD,
	})
	if err != nil {
		return err
	}

	enzymeNames, err := readEnzymeNames(cfg.enzFlag)
	if err != nil {
		return err
	}
	enzymes, err := lookupEnzymes(enzymeNames)
	if err != nil {
		return err
	}

	refBases := design.GenomeBases{}
	genomeBases := cfg.genomeBases
	if genomeBases <= 0 {
		refBases, err = design.CountReferenceBases(cfg.fastaPath)
		if err != nil {
			return err
		}
		switch cfg.denominator {
		case "non-n":
			genomeBases = refBases.NonNBases
		case "all":
			genomeBases = refBases.AllBases
		default:
			return usageError{err: fmt.Errorf("invalid --denominator %q; use non-n or all", cfg.denominator)}
		}
	} else {
		refBases = design.GenomeBases{AllBases: genomeBases, NonNBases: genomeBases}
	}
	if genomeBases <= 0 {
		return fmt.Errorf("genome denominator must be > 0")
	}
	if _, err := fmt.Fprintf(stderr, "genome_bases\t%d\tdenominator\t%s\n", genomeBases, cfg.denominator); err != nil {
		return err
	}
	if err := writeDesignSizeSelectionSummary(stderr, selector.Config()); err != nil {
		return err
	}

	buildWorkers := resolveBuildWorkers(cfg.buildWorkers, cfg.jobs, cfg.threads, len(enzymes))
	idx, err := screen.BuildCutIndexFromFASTAParallel(cfg.fastaPath, enzymes, digest.Options{StrictCuts: cfg.strictCuts}, buildWorkers)
	if err != nil {
		return err
	}
	pairs := idx.PairNames()
	if cfg.maxPairs > 0 && cfg.maxPairs < len(pairs) {
		pairs = pairs[:cfg.maxPairs]
	}
	workers := resolveWorkers(cfg.jobs, cfg.threads, len(pairs))
	if _, err := fmt.Fprintf(stderr, "candidate_enzymes\t%d\n", len(enzymeNames)); err != nil {
		return err
	}
	if _, err := fmt.Fprintf(stderr, "candidate_pairs\t%d\n", len(pairs)); err != nil {
		return err
	}
	if _, err := fmt.Fprintf(stderr, "build_workers\t%d\n", buildWorkers); err != nil {
		return err
	}
	if _, err := fmt.Fprintf(stderr, "pair_score_workers\t%d\n", workers); err != nil {
		return err
	}

	opt := digest.Options{AllowSame: cfg.allowSame, IncludeEnds: cfg.includeEnds, StrictCuts: cfg.strictCuts}
	summaries, err := scorePairs(idx, pairs, selector, opt, workers)
	if err != nil {
		return err
	}

	budget := design.SequencingBudget{
		ReadLayout:           cfg.readLayout,
		ReadLength:           cfg.readLength,
		LaneReadPairs:        cfg.laneReadPairs,
		Lanes:                cfg.lanes,
		UsableReadFraction:   cfg.usableReadFraction,
		Samples:              cfg.samples,
		TargetMeanLocusDepth: cfg.desiredDepth,
	}
	target := design.DesignTarget{
		TargetGenomePct:      cfg.targetGenomePct,
		CoverageTolerancePct: cfg.coverageTolerancePct,
		Objective:            objective,
	}
	weights := design.ScoreWeights{
		Coverage:     cfg.weightCoverage,
		Depth:        cfg.weightDepth,
		Overcoverage: cfg.weightOvercoverage,
		Insert:       cfg.weightInsert,
	}

	candidates := make([]design.Candidate, 0, len(summaries))
	for _, summary := range summaries {
		candidates = append(candidates, design.EvaluateSummary(summary, genomeBases, budget, target, weights))
	}
	design.SortCandidates(candidates, objective)
	reported := candidates
	if cfg.top > 0 && cfg.top < len(reported) {
		reported = append([]design.Candidate(nil), reported[:cfg.top]...)
	}

	warnings := make([]string, 0)
	feasiblePairs := 0
	for _, candidate := range candidates {
		if candidate.Feasible {
			feasiblePairs++
		}
	}
	if feasiblePairs == 0 {
		warnings = append(warnings, "no enzyme pair matched both target coverage tolerance and target mean locus depth under the supplied budget")
	}
	if len(candidates) == 0 {
		warnings = append(warnings, "no candidate pairs were scored")
	}

	tsvPath, jsonPath := resolveOutputPaths(cfg)
	if err := ensureOutputPaths(tsvPath, jsonPath, cfg.force); err != nil {
		return err
	}

	report := buildReport(args, cfg, idx, refBases, genomeBases, selector.Config(), budget, target, weights, warnings, candidates, reported, tsvPath, jsonPath)
	if err := writeCandidatesTSV(tsvPath, report.Results); err != nil {
		return err
	}
	if err := writeJSONAtomic(jsonPath, report); err != nil {
		return err
	}

	if _, err := fmt.Fprintf(stderr, "design_tsv\t%s\n", tsvPath); err != nil {
		return err
	}
	if _, err := fmt.Fprintf(stderr, "design_json\t%s\n", jsonPath); err != nil {
		return err
	}
	return nil
}

func parseArgs(args []string, stderr io.Writer) (cliConfig, error) {
	cfg := cliConfig{}
	defaults := design.DefaultScoreWeights()

	fs := flag.NewFlagSet("radigest-design", flag.ContinueOnError)
	fs.SetOutput(stderr)

	fs.StringVar(&cfg.fastaPath, "fasta", "", "reference FASTA file")
	fs.StringVar(&cfg.fastaPath, "ref", "", "alias for --fasta")
	fs.StringVar(&cfg.enzFlag, "enzymes", "", "comma-separated enzymes, a file with enzyme names, or 'all'")
	fs.StringVar(&cfg.outDir, "out-dir", "radigest_design", "output directory for design.tsv and design.json")
	fs.StringVar(&cfg.tsvPath, "tsv", "", "explicit output TSV path; default <out-dir>/design.tsv")
	fs.StringVar(&cfg.jsonPath, "json", "", "explicit output JSON path; default <out-dir>/design.json")
	fs.StringVar(&cfg.denominator, "denominator", "non-n", "FASTA denominator for genome percentages: non-n or all")
	genomeBasesFlag := fs.String("genome-bases", "", "explicit genome denominator, e.g. 2643888753")

	fs.IntVar(&cfg.minLen, "min", 300, "minimum fragment length (bp) for hard size selection")
	fs.IntVar(&cfg.maxLen, "max", 600, "maximum fragment length (bp) for hard size selection")
	fs.IntVar(&cfg.scoreMin, "score-min", 1, "minimum fragment length included in size-selection scoring")
	fs.IntVar(&cfg.scoreMax, "score-max", 2000, "maximum fragment length included in size-selection scoring")
	fs.StringVar(&cfg.sizeModel, "size-model", "normal", "size-selection model: hard, normal, triangular, or soft-window")
	fs.Float64Var(&cfg.sizeMean, "size-mean", 275, "target/peak insert length for normal/triangular models")
	fs.Float64Var(&cfg.sizeSD, "size-sd", 85, "standard deviation for --size-model normal")
	fs.Float64Var(&cfg.sizeEdgeSD, "size-edge-sd", 25, "edge softness for --size-model soft-window")
	fs.BoolVar(&cfg.allowSame, "allow-same", false, "double digest: also keep AA/BB neighbors (default AB/BA only)")
	fs.BoolVar(&cfg.includeEnds, "include-ends", false, "also score terminal fragments from contig ends to nearest cut")
	fs.BoolVar(&cfg.strictCuts, "strict-cuts", false, "error if an enzyme lacks a caret and CutIndex==0")

	fs.StringVar(&cfg.readLayout, "read-layout", "pe", "sequencing layout for insert diagnostics: pe or se")
	fs.IntVar(&cfg.readLength, "read-length", 0, "read length in bp, e.g. 150")
	laneReadPairsFlag := fs.String("lane-read-pairs", "", "PE read pairs per lane, e.g. 300M")
	flowcellReadPairsFlag := fs.String("flowcell-read-pairs", "", "total read pairs across a one-flowcell/run budget, e.g. 50M")
	fs.IntVar(&cfg.lanes, "lanes", 1, "number of lanes in the sequencing budget")
	fs.Float64Var(&cfg.usableReadFraction, "usable-read-fraction", 1.0, "fraction of reads usable after demultiplexing/QC/deduplication")
	fs.IntVar(&cfg.samples, "samples", 0, "planned number of samples")
	fs.Float64Var(&cfg.desiredDepth, "desired-depth", 0, "target mean read-pair depth per recovered locus")
	fs.Float64Var(&cfg.desiredDepth, "target-depth", 0, "alias for --desired-depth")
	fs.Float64Var(&cfg.desiredDepth, "depth", 0, "alias for --desired-depth")
	fs.Float64Var(&cfg.targetGenomePct, "target-genome-pct", 0, "target weighted genome percentage")
	fs.Float64Var(&cfg.targetGenomePct, "pct", 0, "alias for --target-genome-pct")
	fs.Float64Var(&cfg.coverageTolerancePct, "coverage-tolerance-pct", 0.25, "absolute genome-percentage tolerance around --target-genome-pct")
	fs.StringVar(&cfg.objective, "objective", string(design.ObjectiveBalanced), "ranking objective: balanced, closest-coverage, depth-first, feasible-lowest-coverage, or max-depth")
	fs.Float64Var(&cfg.weightCoverage, "weight-coverage", defaults.Coverage, "fit-loss weight for coverage error")
	fs.Float64Var(&cfg.weightDepth, "weight-depth", defaults.Depth, "fit-loss weight for depth shortfall")
	fs.Float64Var(&cfg.weightOvercoverage, "weight-overcoverage", defaults.Overcoverage, "additional fit-loss weight for overcoverage")
	fs.Float64Var(&cfg.weightInsert, "weight-insert", defaults.Insert, "fit-loss weight for insert-size risk")

	fs.IntVar(&cfg.jobs, "jobs", 0, "parallel pair-scoring workers (default: --threads)")
	fs.IntVar(&cfg.threads, "threads", runtime.NumCPU(), "worker count alias used when --jobs is not set")
	fs.IntVar(&cfg.buildWorkers, "build-workers", 0, "parallel cut-index build workers (default: --jobs, then --threads)")
	fs.IntVar(&cfg.maxPairs, "max-pairs", 0, "maximum number of pairs to score; 0 means all pairs")
	fs.IntVar(&cfg.top, "top", 0, "limit reported rows/results after ranking; 0 means all")
	fs.BoolVar(&cfg.force, "force", false, "overwrite existing TSV/JSON outputs")
	fs.BoolVar(&cfg.showVersion, "version", false, "print version and exit")

	fs.Usage = func() {
		_, _ = fmt.Fprintln(stderr, "radigest-design — inverse enzyme-pair design from coverage/depth targets")
		_, _ = fmt.Fprintln(stderr)
		_, _ = fmt.Fprintln(stderr, "Usage:")
		_, _ = fmt.Fprintln(stderr, "  radigest-design --ref ref.fa --enzymes enzymes.txt --pct 2.5 --depth 10 --samples 96 --read-length 150 --flowcell-read-pairs 300M [options]")
		_, _ = fmt.Fprintln(stderr, "  radigest-design --fasta ref.fa --enzymes enzymes.txt --target-genome-pct 2.5 --desired-depth 10 --samples 96 --read-length 150 --lane-read-pairs 300M [options]")
		_, _ = fmt.Fprintln(stderr)
		fs.PrintDefaults()
	}

	if err := fs.Parse(args); err != nil {
		if errors.Is(err, flag.ErrHelp) {
			cfg.help = true
			return cfg, nil
		}
		return cfg, usageError{err: err}
	}
	lanesExplicit := false
	fs.Visit(func(f *flag.Flag) {
		if f.Name == "lanes" {
			lanesExplicit = true
		}
	})
	if cfg.showVersion {
		return cfg, nil
	}
	if cfg.fastaPath == "" {
		return cfg, usageError{err: errors.New("--fasta/--ref is required")}
	}
	if cfg.enzFlag == "" {
		return cfg, usageError{err: errors.New("--enzymes is required")}
	}
	if cfg.readLength <= 0 {
		return cfg, usageError{err: errors.New("--read-length is required and must be > 0")}
	}
	if cfg.samples <= 0 {
		return cfg, usageError{err: errors.New("--samples is required and must be > 0")}
	}
	if cfg.desiredDepth <= 0 || math.IsNaN(cfg.desiredDepth) || math.IsInf(cfg.desiredDepth, 0) {
		return cfg, usageError{err: errors.New("--desired-depth/--depth is required and must be a finite value > 0")}
	}
	if cfg.targetGenomePct <= 0 || math.IsNaN(cfg.targetGenomePct) || math.IsInf(cfg.targetGenomePct, 0) {
		return cfg, usageError{err: errors.New("--target-genome-pct/--pct is required and must be a finite value > 0")}
	}
	if cfg.coverageTolerancePct < 0 || math.IsNaN(cfg.coverageTolerancePct) || math.IsInf(cfg.coverageTolerancePct, 0) {
		return cfg, usageError{err: fmt.Errorf("--coverage-tolerance-pct must be >= 0 (got %g)", cfg.coverageTolerancePct)}
	}
	if cfg.minLen < 0 || cfg.maxLen < cfg.minLen {
		return cfg, usageError{err: fmt.Errorf("invalid hard size window: min=%d max=%d", cfg.minLen, cfg.maxLen)}
	}
	if cfg.scoreMin < 0 || cfg.scoreMax < cfg.scoreMin {
		return cfg, usageError{err: fmt.Errorf("invalid score window: score-min=%d score-max=%d", cfg.scoreMin, cfg.scoreMax)}
	}
	if cfg.lanes <= 0 {
		return cfg, usageError{err: fmt.Errorf("--lanes must be > 0 (got %d)", cfg.lanes)}
	}
	if cfg.usableReadFraction <= 0 || cfg.usableReadFraction > 1 || math.IsNaN(cfg.usableReadFraction) || math.IsInf(cfg.usableReadFraction, 0) {
		return cfg, usageError{err: fmt.Errorf("--usable-read-fraction must be in (0,1] (got %g)", cfg.usableReadFraction)}
	}
	cfg.readLayout = strings.ToLower(strings.TrimSpace(cfg.readLayout))
	if cfg.readLayout != "pe" && cfg.readLayout != "se" {
		return cfg, usageError{err: fmt.Errorf("invalid --read-layout %q; use pe or se", cfg.readLayout)}
	}
	if *laneReadPairsFlag == "" && *flowcellReadPairsFlag == "" {
		return cfg, usageError{err: errors.New("exactly one of --lane-read-pairs or --flowcell-read-pairs is required")}
	}
	if *laneReadPairsFlag != "" && *flowcellReadPairsFlag != "" {
		return cfg, usageError{err: errors.New("use only one of --lane-read-pairs or --flowcell-read-pairs")}
	}
	if *laneReadPairsFlag != "" {
		value, err := parsePositiveCount(*laneReadPairsFlag)
		if err != nil {
			return cfg, usageError{err: fmt.Errorf("--lane-read-pairs: %w", err)}
		}
		cfg.laneReadPairs = value
	} else {
		if lanesExplicit {
			return cfg, usageError{err: errors.New("--lanes is only meaningful with --lane-read-pairs; omit --lanes when using --flowcell-read-pairs")}
		}
		value, err := parsePositiveCount(*flowcellReadPairsFlag)
		if err != nil {
			return cfg, usageError{err: fmt.Errorf("--flowcell-read-pairs: %w", err)}
		}
		cfg.flowcellReadPairs = value
		cfg.laneReadPairs = value
		cfg.lanes = 1
	}
	if *genomeBasesFlag != "" {
		value, err := parsePositiveCountInt(*genomeBasesFlag)
		if err != nil {
			return cfg, usageError{err: fmt.Errorf("--genome-bases: %w", err)}
		}
		cfg.genomeBases = value
	}
	if cfg.buildWorkers < 0 || cfg.jobs < 0 || cfg.threads < 0 || cfg.maxPairs < 0 || cfg.top < 0 {
		return cfg, usageError{err: errors.New("--build-workers, --jobs, --threads, --max-pairs, and --top must be >= 0")}
	}
	for name, weight := range map[string]float64{
		"--weight-coverage":     cfg.weightCoverage,
		"--weight-depth":        cfg.weightDepth,
		"--weight-overcoverage": cfg.weightOvercoverage,
		"--weight-insert":       cfg.weightInsert,
	} {
		if weight < 0 || math.IsNaN(weight) || math.IsInf(weight, 0) {
			return cfg, usageError{err: fmt.Errorf("%s must be finite and >= 0 (got %g)", name, weight)}
		}
	}
	return cfg, nil
}

func writeDesignSizeSelectionSummary(stderr io.Writer, cfg sizeselect.Config) error {
	if _, err := fmt.Fprintf(stderr, "size_model\t%s\n", cfg.Model); err != nil {
		return err
	}
	if _, err := fmt.Fprintf(stderr, "hard_size_window_bp\t%d-%d\n", cfg.Min, cfg.Max); err != nil {
		return err
	}
	if _, err := fmt.Fprintf(stderr, "score_range_bp\t%d-%d\n", cfg.ScoreMin, cfg.ScoreMax); err != nil {
		return err
	}

	mean := "NA"
	sd := "NA"
	switch cfg.Model {
	case sizeselect.ModelNormal:
		mean = formatDesignStderrFloat(cfg.Mean)
		sd = formatDesignStderrFloat(cfg.SD)
	case sizeselect.ModelTriangular:
		mean = formatDesignStderrFloat(cfg.Mean)
	}
	if _, err := fmt.Fprintf(stderr, "size_mean_bp\t%s\n", mean); err != nil {
		return err
	}
	if _, err := fmt.Fprintf(stderr, "size_sd_bp\t%s\n", sd); err != nil {
		return err
	}
	if cfg.Model == sizeselect.ModelSoftWindow {
		if _, err := fmt.Fprintf(stderr, "size_edge_sd_bp\t%s\n", formatDesignStderrFloat(cfg.EdgeSD)); err != nil {
			return err
		}
	}
	return nil
}

func formatDesignStderrFloat(value float64) string {
	if math.IsNaN(value) || math.IsInf(value, 0) {
		return "NA"
	}
	return strconv.FormatFloat(value, 'f', -1, 64)
}

func readEnzymeNames(value string) ([]string, error) {
	if strings.EqualFold(strings.TrimSpace(value), "all") {
		names := make([]string, 0, len(enzyme.DB))
		for name := range enzyme.DB {
			names = append(names, name)
		}
		sort.Strings(names)
		return names, nil
	}

	var raw []string
	if data, err := os.ReadFile(value); err == nil {
		for _, line := range strings.Split(string(data), "\n") {
			line = strings.TrimSpace(line)
			if line == "" || strings.HasPrefix(line, "#") {
				continue
			}
			raw = append(raw, splitNames(line)...)
		}
	} else {
		raw = splitNames(value)
	}
	deduped := make([]string, 0, len(raw))
	seen := make(map[string]struct{}, len(raw))
	for _, name := range raw {
		name = strings.TrimSpace(name)
		if name == "" {
			continue
		}
		if _, ok := seen[name]; ok {
			continue
		}
		seen[name] = struct{}{}
		deduped = append(deduped, name)
	}
	if len(deduped) < 2 {
		return nil, fmt.Errorf("need at least two candidate enzymes")
	}
	return deduped, nil
}

func splitNames(value string) []string {
	fields := strings.FieldsFunc(value, func(r rune) bool {
		return r == ',' || r == '\t' || r == ' ' || r == '\r'
	})
	out := make([]string, 0, len(fields))
	for _, field := range fields {
		field = strings.TrimSpace(field)
		if field != "" {
			out = append(out, field)
		}
	}
	return out
}

func lookupEnzymes(names []string) ([]enzyme.Enzyme, error) {
	out := make([]enzyme.Enzyme, 0, len(names))
	for _, name := range names {
		enz, ok := enzyme.Get(name)
		if !ok {
			return nil, fmt.Errorf("unknown enzyme %q", name)
		}
		out = append(out, enz)
	}
	return out, nil
}

func resolveWorkers(jobs, threads, pairCount int) int {
	workers := jobs
	if workers <= 0 {
		workers = threads
	}
	if workers <= 0 {
		workers = 1
	}
	if pairCount > 0 && workers > pairCount {
		workers = pairCount
	}
	if workers <= 0 {
		workers = 1
	}
	return workers
}

func resolveBuildWorkers(buildWorkers, jobs, threads, enzymeCount int) int {
	workers := buildWorkers
	if workers <= 0 {
		workers = jobs
	}
	if workers <= 0 {
		workers = threads
	}
	if workers <= 0 {
		workers = 1
	}
	if enzymeCount > 0 && workers > enzymeCount {
		workers = enzymeCount
	}
	if workers <= 0 {
		workers = 1
	}
	return workers
}

func scorePairs(idx screen.CutIndex, pairs []screen.Pair, selector sizeselect.Selector, opt digest.Options, workers int) ([]screen.PairSummary, error) {
	if len(pairs) == 0 {
		return nil, nil
	}
	if workers < 1 {
		workers = 1
	}
	type job struct {
		idx  int
		pair screen.Pair
	}
	type result struct {
		idx     int
		summary screen.PairSummary
		err     error
	}
	jobCh := make(chan job)
	resultCh := make(chan result, len(pairs))
	var wg sync.WaitGroup
	for w := 0; w < workers; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for j := range jobCh {
				summary, err := screen.ScorePair(idx, j.pair.A, j.pair.B, selector, opt)
				resultCh <- result{idx: j.idx, summary: summary, err: err}
			}
		}()
	}
	for i, pair := range pairs {
		jobCh <- job{idx: i, pair: pair}
	}
	close(jobCh)
	wg.Wait()
	close(resultCh)

	summaries := make([]screen.PairSummary, len(pairs))
	var firstErr error
	for res := range resultCh {
		if res.err != nil && firstErr == nil {
			firstErr = res.err
		}
		summaries[res.idx] = res.summary
	}
	if firstErr != nil {
		return nil, firstErr
	}
	return summaries, nil
}

func buildReport(args []string, cfg cliConfig, idx screen.CutIndex, refBases design.GenomeBases, genomeBases int64, selectorCfg sizeselect.Config, budget design.SequencingBudget, target design.DesignTarget, weights design.ScoreWeights, warnings []string, allCandidates []design.Candidate, reported []design.Candidate, tsvPath, jsonPath string) designReport {
	feasiblePairs := 0
	for _, candidate := range allCandidates {
		if candidate.Feasible {
			feasiblePairs++
		}
	}
	summary := runSummary{
		CandidateEnzymes: len(idx.EnzymeNames),
		CandidatePairs:   len(allCandidates),
		ReportedPairs:    len(reported),
		FeasiblePairs:    feasiblePairs,
	}
	if len(allCandidates) > 0 {
		summary.BestPair = []string{allCandidates[0].EnzymeA, allCandidates[0].EnzymeB}
	}
	command := make([]string, 0, len(args)+1)
	command = append(command, "radigest-design")
	command = append(command, args...)

	digestParams := digestParameters{
		MinLength:   cfg.minLen,
		MaxLength:   cfg.maxLen,
		ScoreMin:    selectorCfg.ScoreMin,
		ScoreMax:    selectorCfg.ScoreMax,
		SizeModel:   string(selectorCfg.Model),
		AllowSame:   cfg.allowSame,
		IncludeEnds: cfg.includeEnds,
		StrictCuts:  cfg.strictCuts,
	}
	switch selectorCfg.Model {
	case sizeselect.ModelNormal:
		digestParams.SizeMean = selectorCfg.Mean
		digestParams.SizeSD = selectorCfg.SD
	case sizeselect.ModelTriangular:
		digestParams.SizeMean = selectorCfg.Mean
	case sizeselect.ModelSoftWindow:
		digestParams.SizeEdgeSD = selectorCfg.EdgeSD
	}

	return designReport{
		SchemaVersion:   design.SchemaVersion,
		RadigestVersion: version,
		Command:         command,
		Input: inputSummary{
			FASTA:       cfg.fastaPath,
			Denominator: cfg.denominator,
			GenomeBases: genomeBases,
			Reference:   refBases,
		},
		Digest:     digestParams,
		Sequencing: budget,
		Target:     target,
		Weights:    weights,
		Screening: screen.ScreeningStats{
			Engine:                   screen.EngineCachedCutIndex,
			CandidateEnzymes:         len(idx.EnzymeNames),
			Records:                  len(idx.Records),
			CachedCutSites:           idx.CachedCutSites(),
			CacheMemoryEstimateBytes: idx.CacheMemoryEstimateBytes(),
		},
		Outputs:  outputSummary{TSV: tsvPath, JSON: jsonPath},
		Warnings: warnings,
		Summary:  summary,
		Results:  append([]design.Candidate(nil), reported...),
	}
}

func resolveOutputPaths(cfg cliConfig) (string, string) {
	tsvPath := strings.TrimSpace(cfg.tsvPath)
	jsonPath := strings.TrimSpace(cfg.jsonPath)
	if tsvPath == "" {
		tsvPath = filepath.Join(cfg.outDir, "design.tsv")
	}
	if jsonPath == "" {
		jsonPath = filepath.Join(cfg.outDir, "design.json")
	}
	return tsvPath, jsonPath
}

func ensureOutputPaths(tsvPath, jsonPath string, force bool) error {
	for _, path := range []string{tsvPath, jsonPath} {
		if strings.TrimSpace(path) == "" {
			return fmt.Errorf("output path must not be empty")
		}
		if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
			return err
		}
		if !force {
			if _, err := os.Stat(path); err == nil {
				return fmt.Errorf("output exists: %s; use --force to overwrite", path)
			} else if !os.IsNotExist(err) {
				return err
			}
		}
	}
	if samePath(tsvPath, jsonPath) {
		return fmt.Errorf("refusing to write TSV and JSON to the same path %q", tsvPath)
	}
	return nil
}

func samePath(a, b string) bool {
	absA, errA := filepath.Abs(a)
	absB, errB := filepath.Abs(b)
	if errA != nil || errB != nil {
		return a == b
	}
	return filepath.Clean(absA) == filepath.Clean(absB)
}

func writeCandidatesTSV(path string, candidates []design.Candidate) (err error) {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer func() {
		if closeErr := f.Close(); closeErr != nil && err == nil {
			err = closeErr
		}
	}()
	writer := csv.NewWriter(f)
	writer.Comma = '\t'
	if err := writer.Write(designTSVHeader()); err != nil {
		return err
	}
	for _, candidate := range candidates {
		if err := writer.Write(candidateTSVRow(candidate)); err != nil {
			return err
		}
	}
	writer.Flush()
	return writer.Error()
}

func designTSVHeader() []string {
	return []string{
		"rank",
		"enzyme_a",
		"enzyme_b",
		"feasible",
		"decision_reason",
		"fit_score",
		"fit_loss",
		"target_genome_pct",
		"predicted_weighted_genome_pct",
		"coverage_error_pct_points",
		"coverage_error_rel",
		"overcoverage_rel",
		"undercoverage_rel",
		"target_mean_locus_depth",
		"predicted_mean_locus_depth",
		"depth_margin",
		"depth_shortfall_rel",
		"read_pairs_per_sample",
		"required_pairs_per_sample_full_target",
		"weighted_bases",
		"weighted_fragments",
		"mean_weighted_length",
		"raw_bases_in_window",
		"raw_fragments_in_window",
		"budget_supported_genome_pct",
		"budget_supported_weighted_bases",
		"max_samples_per_lane_full_target",
		"max_samples_total_full_target",
		"lanes_required_full_target",
		"adapter_threshold_bp",
		"overlap_threshold_bp",
		"mean_insert_category",
		"insert_penalty",
		"records",
		"cached_cut_sites",
		"cache_memory_estimate_bytes",
	}
}

func candidateTSVRow(c design.Candidate) []string {
	return []string{
		strconv.Itoa(c.Rank),
		c.EnzymeA,
		c.EnzymeB,
		strconv.FormatBool(c.Feasible),
		c.DecisionReason,
		formatFloat(c.FitScore),
		formatFloat(c.FitLoss),
		formatFloat(c.TargetGenomePct),
		formatFloat(c.PredictedWeightedGenomePct),
		formatFloat(c.CoverageErrorPctPoints),
		formatFloat(c.CoverageErrorRel),
		formatFloat(c.OvercoverageRel),
		formatFloat(c.UndercoverageRel),
		formatFloat(c.TargetMeanLocusDepth),
		formatFloat(c.PredictedMeanLocusDepth),
		formatFloat(c.DepthMargin),
		formatFloat(c.DepthShortfallRel),
		formatFloat(c.ReadPairsPerSample),
		formatFloat(c.RequiredPairsPerSampleTarget),
		formatFloat(c.WeightedBases),
		formatFloat(c.WeightedFragments),
		formatFloat(c.MeanWeightedLength),
		strconv.FormatInt(c.RawBasesInWindow, 10),
		strconv.Itoa(c.RawFragmentsInWindow),
		formatFloat(c.BudgetSupportedGenomePct),
		formatFloat(c.BudgetSupportedWeightedBases),
		strconv.Itoa(c.MaxSamplesPerLaneFullTarget),
		strconv.Itoa(c.MaxSamplesTotalFullTarget),
		strconv.Itoa(c.LanesRequiredFullTarget),
		strconv.Itoa(c.AdapterThresholdBP),
		zeroIntBlank(c.OverlapThresholdBP),
		c.MeanInsertCategory,
		formatFloat(c.InsertPenalty),
		strconv.Itoa(c.Records),
		strconv.Itoa(c.CachedCutSites),
		strconv.FormatInt(c.CacheMemoryEstimateBytes, 10),
	}
}

func writeJSONAtomic(path string, value any) error {
	dir := filepath.Dir(path)
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return err
	}
	tmp, err := os.CreateTemp(dir, filepath.Base(path)+".*.tmp")
	if err != nil {
		return err
	}
	tmpName := tmp.Name()
	keepTmp := false
	defer func() {
		if !keepTmp {
			_ = os.Remove(tmpName)
		}
	}()
	if err := tmp.Chmod(0o644); err != nil {
		_ = tmp.Close()
		return err
	}
	enc := json.NewEncoder(tmp)
	enc.SetIndent("", "  ")
	if err := enc.Encode(value); err != nil {
		_ = tmp.Close()
		return err
	}
	if err := tmp.Sync(); err != nil {
		_ = tmp.Close()
		return err
	}
	if err := tmp.Close(); err != nil {
		return err
	}
	if err := os.Rename(tmpName, path); err != nil {
		return err
	}
	keepTmp = true
	return nil
}

func parsePositiveCount(value string) (float64, error) {
	text := strings.TrimSpace(strings.ReplaceAll(value, "_", ""))
	if text == "" {
		return 0, fmt.Errorf("must be a finite count > 0")
	}
	multiplier := 1.0
	suffix := strings.ToLower(text[len(text)-1:])
	switch suffix {
	case "k":
		multiplier = 1e3
		text = text[:len(text)-1]
	case "m":
		multiplier = 1e6
		text = text[:len(text)-1]
	case "g":
		multiplier = 1e9
		text = text[:len(text)-1]
	case "t":
		multiplier = 1e12
		text = text[:len(text)-1]
	}
	if text == "" {
		return 0, fmt.Errorf("must be a finite count > 0")
	}
	parsed, err := strconv.ParseFloat(text, 64)
	if err != nil {
		return 0, fmt.Errorf("invalid count %q", value)
	}
	parsed *= multiplier
	if parsed <= 0 || math.IsNaN(parsed) || math.IsInf(parsed, 0) {
		return 0, fmt.Errorf("must be a finite count > 0")
	}
	return parsed, nil
}

func parsePositiveCountInt(value string) (int64, error) {
	parsed, err := parsePositiveCount(value)
	if err != nil {
		return 0, err
	}
	rounded := math.Round(parsed)
	if math.Abs(parsed-rounded) > 1e-6 {
		return 0, fmt.Errorf("must resolve to an integer count")
	}
	return int64(rounded), nil
}

func formatFloat(value float64) string {
	if math.IsNaN(value) || math.IsInf(value, 0) {
		return ""
	}
	return strconv.FormatFloat(value, 'f', 6, 64)
}

func zeroIntBlank(value int) string {
	if value == 0 {
		return ""
	}
	return strconv.Itoa(value)
}
