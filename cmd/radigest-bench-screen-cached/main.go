package main

import (
	"bytes"
	"encoding/csv"
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"hash/fnv"
	"io"
	"os"
	"path/filepath"
	"regexp"
	"runtime"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/ericksamera/radigest/internal/digest"
	"github.com/ericksamera/radigest/internal/enzyme"
	"github.com/ericksamera/radigest/internal/screen"
	"github.com/ericksamera/radigest/internal/sizeselect"
)

var version = "dev"

var safeTag = regexp.MustCompile(`[^A-Za-z0-9_.-]+`)

type usageError struct{ err error }

func (e usageError) Error() string { return e.err.Error() }
func (e usageError) Unwrap() error { return e.err }

type outputMode string

const (
	outputModeNone    outputMode = "none"
	outputModeMarshal outputMode = "marshal"
	outputModeWrite   outputMode = "write"
)

type benchConfig struct {
	fastaPath    string
	enzFlag      string
	outDir       string
	minLen       int
	maxLen       int
	scoreMin     int
	scoreMax     int
	sizeModel    string
	sizeMean     float64
	sizeSD       float64
	sizeEdgeSD   float64
	jobs         int
	threads      int
	buildWorkers int
	runs         int
	maxPairs     int
	allowSame    bool
	includeEnds  bool
	strictCuts   bool
	reuseIndex   bool
	mode         outputMode
	indentJSON   bool
	syncJSON     bool
	writeLogs    bool
	force        bool
	showVersion  bool
	help         bool
}

type benchRun struct {
	Run                      int     `json:"run"`
	CandidateEnzymes         int     `json:"candidate_enzymes"`
	CandidatePairs           int     `json:"candidate_pairs"`
	Jobs                     int     `json:"jobs"`
	GOMAXPROCS               int     `json:"gomaxprocs"`
	BuildWorkers             int     `json:"build_workers"`
	Records                  int     `json:"records"`
	CachedCutSites           int     `json:"cached_cut_sites"`
	CacheMemoryEstimateBytes int64   `json:"cache_memory_estimate_bytes"`
	BuildCutIndexSeconds     float64 `json:"build_cut_index_seconds"`
	ScorePairsSeconds        float64 `json:"score_pairs_seconds"`
	JSONMarshalSeconds       float64 `json:"json_marshal_seconds"`
	WriteJSONSeconds         float64 `json:"write_json_seconds"`
	TotalSeconds             float64 `json:"total_seconds"`
	PairsPerSecondScorePhase float64 `json:"pairs_per_second_score_phase"`
	PairsPerSecondEndToEnd   float64 `json:"pairs_per_second_end_to_end"`
	Summaries                int     `json:"summaries"`
	TotalFragments           int64   `json:"total_fragments"`
	TotalBases               int64   `json:"total_bases"`
	DigestChecksum           string  `json:"digest_checksum"`
	OutputBytes              int64   `json:"output_bytes"`
	OutputMode               string  `json:"output_mode"`
	SyncJSON                 bool    `json:"sync_json"`
	WriteLogs                bool    `json:"write_logs"`
	ReusedIndex              bool    `json:"reused_index"`
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
		_, err := fmt.Fprintf(stdout, "radigest-bench-screen-cached %s\n", version)
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
	opt := digest.Options{
		AllowSame:   cfg.allowSame,
		IncludeEnds: cfg.includeEnds,
		StrictCuts:  cfg.strictCuts,
	}

	buildWorkers := resolveBuildWorkers(cfg.buildWorkers, cfg.jobs, cfg.threads, len(enzymes))

	var reusedIndex screen.CutIndex
	var reusedBuildSeconds float64
	if cfg.reuseIndex {
		start := time.Now()
		idx, err := screen.BuildCutIndexFromFASTAParallel(cfg.fastaPath, enzymes, digest.Options{StrictCuts: cfg.strictCuts}, buildWorkers)
		if err != nil {
			return err
		}
		reusedIndex = idx
		reusedBuildSeconds = time.Since(start).Seconds()
		if _, err := fmt.Fprintf(stderr, "prebuilt_cut_index_seconds\t%.6f\n", reusedBuildSeconds); err != nil {
			return err
		}
	}

	writer := csv.NewWriter(stdout)
	writer.Comma = '\t'
	if err := writer.Write(benchHeader()); err != nil {
		return err
	}

	for runID := 1; runID <= cfg.runs; runID++ {
		result, err := runBenchmarkOnce(runID, cfg, enzymes, selector, opt, reusedIndex, cfg.reuseIndex, reusedBuildSeconds, buildWorkers)
		if err != nil {
			return err
		}
		if err := writer.Write(result.tsvRow()); err != nil {
			return err
		}
		writer.Flush()
		if err := writer.Error(); err != nil {
			return err
		}
	}
	return nil
}

func parseArgs(args []string, stderr io.Writer) (benchConfig, error) {
	cfg := benchConfig{}

	fs := flag.NewFlagSet("radigest-bench-screen-cached", flag.ContinueOnError)
	fs.SetOutput(stderr)

	fs.StringVar(&cfg.fastaPath, "fasta", "", "reference FASTA file")
	fs.StringVar(&cfg.enzFlag, "enzymes", "", "comma-separated enzymes or a file with one/comma-separated enzymes per line")
	fs.StringVar(&cfg.outDir, "out-dir", "pair_screen_bench", "output directory used only with --output-mode write")
	fs.IntVar(&cfg.minLen, "min", 300, "minimum fragment length (bp) for hard size selection")
	fs.IntVar(&cfg.maxLen, "max", 600, "maximum fragment length (bp) for hard size selection")
	fs.IntVar(&cfg.scoreMin, "score-min", 1, "minimum fragment length included in size-selection scoring")
	fs.IntVar(&cfg.scoreMax, "score-max", 2000, "maximum fragment length included in size-selection scoring")
	fs.StringVar(&cfg.sizeModel, "size-model", "normal", "size-selection model: hard, normal, triangular, or soft-window")
	fs.Float64Var(&cfg.sizeMean, "size-mean", 275, "target/peak insert length for normal/triangular models")
	fs.Float64Var(&cfg.sizeSD, "size-sd", 85, "standard deviation for --size-model normal")
	fs.Float64Var(&cfg.sizeEdgeSD, "size-edge-sd", 25, "edge softness for --size-model soft-window")
	fs.IntVar(&cfg.jobs, "jobs", 0, "parallel pair-scoring workers (default: --threads)")
	fs.IntVar(&cfg.threads, "threads", runtime.NumCPU(), "worker count alias used when --jobs is not set")
	fs.IntVar(&cfg.buildWorkers, "build-workers", 0, "parallel cut-index build workers (default: --jobs, then --threads); scans candidate enzymes concurrently per FASTA record")
	fs.IntVar(&cfg.runs, "runs", 1, "number of benchmark repetitions")
	fs.IntVar(&cfg.maxPairs, "max-pairs", 0, "maximum number of pairs to score; 0 means all pairs")
	fs.BoolVar(&cfg.allowSame, "allow-same", false, "double digest: also keep AA/BB neighbors (default AB/BA only)")
	fs.BoolVar(&cfg.includeEnds, "include-ends", false, "also score terminal fragments from contig ends to nearest cut")
	fs.BoolVar(&cfg.strictCuts, "strict-cuts", false, "error if an enzyme lacks a caret and CutIndex==0")
	fs.BoolVar(&cfg.reuseIndex, "reuse-index", false, "build the cut index once and reuse it across --runs")
	mode := fs.String("output-mode", string(outputModeNone), "output phase to include: none, marshal, or write")
	fs.BoolVar(&cfg.indentJSON, "indent-json", false, "use indented JSON in marshal/write output phases")
	fs.BoolVar(&cfg.syncJSON, "sync-json", true, "write JSON atomically and fsync the temporary file in --output-mode write")
	fs.BoolVar(&cfg.writeLogs, "write-logs", true, "write one per-pair log file in --output-mode write")
	fs.BoolVar(&cfg.force, "force", false, "overwrite existing JSON files in --output-mode write")
	fs.BoolVar(&cfg.showVersion, "version", false, "print version and exit")

	fs.Usage = func() {
		_, _ = fmt.Fprintln(stderr, "radigest-bench-screen-cached — phase-timed cached pair-screen benchmark")
		_, _ = fmt.Fprintln(stderr)
		_, _ = fmt.Fprintln(stderr, "Usage:")
		_, _ = fmt.Fprintln(stderr, "  radigest-bench-screen-cached --fasta ref.fa --enzymes enzymes.txt --jobs 4 [options]")
		_, _ = fmt.Fprintln(stderr)
		_, _ = fmt.Fprintln(stderr, "Output is TSV on stdout with one row per run.")
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
	if cfg.showVersion {
		return cfg, nil
	}
	if cfg.fastaPath == "" {
		return cfg, usageError{err: errors.New("--fasta is required")}
	}
	if cfg.enzFlag == "" {
		return cfg, usageError{err: errors.New("--enzymes is required")}
	}
	if cfg.minLen < 0 || cfg.maxLen < cfg.minLen {
		return cfg, usageError{err: fmt.Errorf("invalid hard size window: min=%d max=%d", cfg.minLen, cfg.maxLen)}
	}
	if cfg.runs < 1 {
		return cfg, usageError{err: fmt.Errorf("--runs must be >= 1 (got %d)", cfg.runs)}
	}
	if cfg.buildWorkers < 0 {
		return cfg, usageError{err: fmt.Errorf("--build-workers must be >= 0 (got %d)", cfg.buildWorkers)}
	}
	if cfg.maxPairs < 0 {
		return cfg, usageError{err: fmt.Errorf("--max-pairs must be >= 0 (got %d)", cfg.maxPairs)}
	}
	cfg.mode = outputMode(strings.TrimSpace(strings.ToLower(*mode)))
	switch cfg.mode {
	case outputModeNone, outputModeMarshal, outputModeWrite:
		// ok
	default:
		return cfg, usageError{err: fmt.Errorf("invalid --output-mode %q; expected none, marshal, or write", *mode)}
	}
	return cfg, nil
}

func runBenchmarkOnce(runID int, cfg benchConfig, enzymes []enzyme.Enzyme, selector sizeselect.Selector, opt digest.Options, reusedIndex screen.CutIndex, reuseIndex bool, reusedBuildSeconds float64, buildWorkers int) (benchRun, error) {
	totalStart := time.Now()

	idx := reusedIndex
	buildSeconds := reusedBuildSeconds
	if !reuseIndex {
		start := time.Now()
		var err error
		idx, err = screen.BuildCutIndexFromFASTAParallel(cfg.fastaPath, enzymes, digest.Options{StrictCuts: cfg.strictCuts}, buildWorkers)
		if err != nil {
			return benchRun{}, err
		}
		buildSeconds = time.Since(start).Seconds()
	}

	pairs := idx.PairNames()
	if cfg.maxPairs > 0 && cfg.maxPairs < len(pairs) {
		pairs = pairs[:cfg.maxPairs]
	}
	workers := resolveWorkers(cfg.jobs, cfg.threads, len(pairs))

	scoreStart := time.Now()
	summaries, err := scorePairs(idx, pairs, selector, opt, workers)
	if err != nil {
		return benchRun{}, err
	}
	scoreSeconds := time.Since(scoreStart).Seconds()

	outputBytes := int64(0)
	marshalSeconds := 0.0
	writeSeconds := 0.0
	switch cfg.mode {
	case outputModeMarshal:
		start := time.Now()
		bytesWritten, err := marshalSummaries(summaries, cfg.indentJSON)
		if err != nil {
			return benchRun{}, err
		}
		marshalSeconds = time.Since(start).Seconds()
		outputBytes = bytesWritten
	case outputModeWrite:
		start := time.Now()
		bytesWritten, err := writeSummaries(runID, cfg.outDir, summaries, cfg.indentJSON, cfg.syncJSON, cfg.writeLogs, cfg.force)
		if err != nil {
			return benchRun{}, err
		}
		writeSeconds = time.Since(start).Seconds()
		outputBytes = bytesWritten
	}

	result := summarizeRun(runID, cfg, idx, pairs, summaries)
	result.BuildCutIndexSeconds = buildSeconds
	result.ScorePairsSeconds = scoreSeconds
	result.JSONMarshalSeconds = marshalSeconds
	result.WriteJSONSeconds = writeSeconds
	result.TotalSeconds = time.Since(totalStart).Seconds()
	result.Jobs = workers
	result.GOMAXPROCS = runtime.GOMAXPROCS(0)
	result.BuildWorkers = buildWorkers
	result.OutputMode = string(cfg.mode)
	result.SyncJSON = cfg.syncJSON
	result.WriteLogs = cfg.writeLogs
	result.ReusedIndex = reuseIndex
	result.OutputBytes = outputBytes
	if scoreSeconds > 0 {
		result.PairsPerSecondScorePhase = float64(len(pairs)) / scoreSeconds
	}
	if result.TotalSeconds > 0 {
		result.PairsPerSecondEndToEnd = float64(len(pairs)) / result.TotalSeconds
	}
	return result, nil
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

func marshalSummaries(summaries []screen.PairSummary, indent bool) (int64, error) {
	var total int64
	for _, summary := range summaries {
		var data []byte
		var err error
		if indent {
			data, err = json.MarshalIndent(summary, "", "  ")
		} else {
			data, err = json.Marshal(summary)
		}
		if err != nil {
			return 0, err
		}
		total += int64(len(data) + 1)
	}
	return total, nil
}

func writeSummaries(runID int, outDir string, summaries []screen.PairSummary, indent bool, syncJSON bool, writeLogs bool, force bool) (int64, error) {
	runDir := filepath.Join(outDir, fmt.Sprintf("run_%03d", runID))
	jsonDir := filepath.Join(runDir, "json")
	logDir := filepath.Join(runDir, "logs")
	if err := os.MkdirAll(jsonDir, 0o755); err != nil {
		return 0, err
	}
	if writeLogs {
		if err := os.MkdirAll(logDir, 0o755); err != nil {
			return 0, err
		}
	}
	var total int64
	for _, summary := range summaries {
		tag := pairTag(summary.Enzymes[0], summary.Enzymes[1])
		path := filepath.Join(jsonDir, tag+".json")
		if !force {
			if _, err := os.Stat(path); err == nil {
				return total, fmt.Errorf("output exists: %s; use --force to overwrite", path)
			} else if !os.IsNotExist(err) {
				return total, err
			}
		}
		data, err := encodeSummary(summary, indent)
		if err != nil {
			return total, err
		}
		if err := writeJSONFile(path, data, syncJSON); err != nil {
			return total, err
		}
		total += int64(len(data))
		if writeLogs {
			logData := []byte(fmt.Sprintf("pair\t%s,%s\njson\t%s\nstatus\tok\n", summary.Enzymes[0], summary.Enzymes[1], path))
			logPath := filepath.Join(logDir, tag+".log")
			if err := os.WriteFile(logPath, logData, 0o644); err != nil {
				return total, err
			}
			total += int64(len(logData))
		}
	}
	return total, nil
}

func writeJSONFile(path string, data []byte, syncJSON bool) error {
	if !syncJSON {
		return os.WriteFile(path, data, 0o644)
	}
	dir := filepath.Dir(path)
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
	if _, err := tmp.Write(data); err != nil {
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

func encodeSummary(summary screen.PairSummary, indent bool) ([]byte, error) {
	if indent {
		data, err := json.MarshalIndent(summary, "", "  ")
		if err != nil {
			return nil, err
		}
		return append(data, '\n'), nil
	}
	var buf bytes.Buffer
	if err := json.NewEncoder(&buf).Encode(summary); err != nil {
		return nil, err
	}
	return buf.Bytes(), nil
}

func pairTag(enzymeA, enzymeB string) string {
	return safeTag.ReplaceAllString(enzymeA+"__"+enzymeB, "_")
}

func summarizeRun(runID int, cfg benchConfig, idx screen.CutIndex, pairs []screen.Pair, summaries []screen.PairSummary) benchRun {
	h := fnv.New64a()
	var totalFragments int64
	var totalBases int64
	for _, summary := range summaries {
		totalFragments += int64(summary.TotalFragments)
		totalBases += int64(summary.TotalBases)
		_, _ = fmt.Fprintf(h, "%s,%s:%d:%d:%d:%d;", summary.Enzymes[0], summary.Enzymes[1], summary.TotalFragments, summary.TotalBases, summary.SizeSelection.RawFragmentsScored, summary.SizeSelection.RawBasesScored)
	}
	return benchRun{
		Run:                      runID,
		CandidateEnzymes:         len(idx.EnzymeNames),
		CandidatePairs:           len(pairs),
		Records:                  len(idx.Records),
		CachedCutSites:           idx.CachedCutSites(),
		CacheMemoryEstimateBytes: idx.CacheMemoryEstimateBytes(),
		Summaries:                len(summaries),
		TotalFragments:           totalFragments,
		TotalBases:               totalBases,
		DigestChecksum:           fmt.Sprintf("%016x", h.Sum64()),
	}
}

func readEnzymeNames(value string) ([]string, error) {
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

func benchHeader() []string {
	return []string{
		"run",
		"candidate_enzymes",
		"candidate_pairs",
		"jobs",
		"gomaxprocs",
		"build_workers",
		"records",
		"cached_cut_sites",
		"cache_memory_estimate_bytes",
		"build_cut_index_seconds",
		"score_pairs_seconds",
		"json_marshal_seconds",
		"write_json_seconds",
		"total_seconds",
		"pairs_per_second_score_phase",
		"pairs_per_second_end_to_end",
		"summaries",
		"total_fragments",
		"total_bases",
		"digest_checksum",
		"output_bytes",
		"output_mode",
		"sync_json",
		"write_logs",
		"reused_index",
	}
}

func (r benchRun) tsvRow() []string {
	return []string{
		strconv.Itoa(r.Run),
		strconv.Itoa(r.CandidateEnzymes),
		strconv.Itoa(r.CandidatePairs),
		strconv.Itoa(r.Jobs),
		strconv.Itoa(r.GOMAXPROCS),
		strconv.Itoa(r.BuildWorkers),
		strconv.Itoa(r.Records),
		strconv.Itoa(r.CachedCutSites),
		strconv.FormatInt(r.CacheMemoryEstimateBytes, 10),
		formatFloat(r.BuildCutIndexSeconds),
		formatFloat(r.ScorePairsSeconds),
		formatFloat(r.JSONMarshalSeconds),
		formatFloat(r.WriteJSONSeconds),
		formatFloat(r.TotalSeconds),
		formatFloat(r.PairsPerSecondScorePhase),
		formatFloat(r.PairsPerSecondEndToEnd),
		strconv.Itoa(r.Summaries),
		strconv.FormatInt(r.TotalFragments, 10),
		strconv.FormatInt(r.TotalBases, 10),
		r.DigestChecksum,
		strconv.FormatInt(r.OutputBytes, 10),
		r.OutputMode,
		strconv.FormatBool(r.SyncJSON),
		strconv.FormatBool(r.WriteLogs),
		strconv.FormatBool(r.ReusedIndex),
	}
}

func formatFloat(v float64) string {
	return strconv.FormatFloat(v, 'f', 9, 64)
}
