package main

import (
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"regexp"
	"runtime"
	"sort"
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

type pairJob struct {
	enzymeA string
	enzymeB string
	tag     string
	json    string
	log     string
}

type pairResult struct {
	job pairJob
	err error
}

type pairJSON struct {
	SchemaVersion   int                           `json:"schema_version"`
	RadigestVersion string                        `json:"radigest_version,omitempty"`
	Command         []string                      `json:"command,omitempty"`
	Enzymes         []string                      `json:"enzymes"`
	MinLength       int                           `json:"min_length"`
	MaxLength       int                           `json:"max_length"`
	TotalFragments  int                           `json:"total_fragments"`
	TotalBases      int                           `json:"total_bases"`
	PerChromosome   map[string]screen.RecordStats `json:"per_chromosome"`
	SizeSelection   sizeselect.Stats              `json:"size_selection"`
	Screening       screen.ScreeningStats         `json:"screening"`
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

	fs := flag.NewFlagSet("radigest-screen-pairs-cached", flag.ContinueOnError)
	fs.SetOutput(stderr)

	fastaPath := fs.String("fasta", "", "reference FASTA file")
	enzFlag := fs.String("enzymes", "", "comma-separated enzymes or a file with one/comma-separated enzymes per line")
	outDir := fs.String("out-dir", "pair_screen_cached", "output directory containing json/ and logs/")
	minLen := fs.Int("min", 300, "minimum fragment length (bp) for hard size selection")
	maxLen := fs.Int("max", 600, "maximum fragment length (bp) for hard size selection")
	scoreMin := fs.Int("score-min", 1, "minimum fragment length included in size-selection scoring")
	scoreMax := fs.Int("score-max", 2000, "maximum fragment length included in size-selection scoring")
	sizeModel := fs.String("size-model", "normal", "size-selection model: hard, normal, triangular, or soft-window")
	sizeMean := fs.Float64("size-mean", 275, "target/peak insert length for normal/triangular models")
	sizeSD := fs.Float64("size-sd", 85, "standard deviation for -size-model normal")
	sizeEdgeSD := fs.Float64("size-edge-sd", 25, "edge softness for -size-model soft-window")
	jobsFlag := fs.Int("jobs", 0, "parallel pair-scoring workers (default: -threads)")
	threadsFlag := fs.Int("threads", runtime.NumCPU(), "worker count alias used when -jobs is not set")
	allowSame := fs.Bool("allow-same", false, "double digest: also keep AA/BB neighbors (default AB/BA only)")
	includeEnds := fs.Bool("include-ends", false, "also score terminal fragments from contig ends to nearest cut")
	strictCuts := fs.Bool("strict-cuts", false, "error if an enzyme lacks a caret and CutIndex==0")
	force := fs.Bool("force", false, "overwrite existing JSON files")
	dryRun := fs.Bool("dry-run", false, "print pair jobs without reading FASTA or writing JSON")
	maxPairs := fs.Int("max-pairs", 0, "maximum number of pairs to score; 0 means all pairs")
	showVersion := fs.Bool("version", false, "print version and exit")

	fs.Usage = func() {
		_, _ = fmt.Fprintln(stderr, "radigest-screen-pairs-cached — cached cut-index enzyme-pair screen")
		_, _ = fmt.Fprintln(stderr)
		_, _ = fmt.Fprintln(stderr, "Usage:")
		_, _ = fmt.Fprintln(stderr, "  radigest-screen-pairs-cached --fasta ref.fa.gz --enzymes enzymes.txt [options]")
		_, _ = fmt.Fprintln(stderr)
		fs.PrintDefaults()
	}

	if err := fs.Parse(args); err != nil {
		if errors.Is(err, flag.ErrHelp) {
			return nil
		}
		return usageError{err: err}
	}
	if *showVersion {
		_, err := fmt.Fprintf(stdout, "radigest-screen-pairs-cached %s\n", version)
		return err
	}
	if *fastaPath == "" {
		return usageError{err: errors.New("--fasta is required")}
	}
	if *enzFlag == "" {
		return usageError{err: errors.New("--enzymes is required")}
	}
	if *minLen < 0 || *maxLen < *minLen {
		return usageError{err: fmt.Errorf("invalid hard size window: min=%d max=%d", *minLen, *maxLen)}
	}

	enzymeNames, err := readEnzymeNames(*enzFlag)
	if err != nil {
		return err
	}
	enzymes, err := lookupEnzymes(enzymeNames)
	if err != nil {
		return err
	}
	pairJobs := buildPairJobs(enzymeNames, *outDir, *force, *maxPairs)
	if _, err := fmt.Fprintf(stderr, "candidate_enzymes\t%d\n", len(enzymeNames)); err != nil {
		return err
	}
	if _, err := fmt.Fprintf(stderr, "pair_jobs\t%d\n", len(pairJobs)); err != nil {
		return err
	}
	if len(pairJobs) == 0 {
		return nil
	}
	if *dryRun {
		for _, job := range pairJobs {
			if _, err := fmt.Fprintf(stdout, "%s\t%s,%s\t%s\n", job.tag, job.enzymeA, job.enzymeB, job.json); err != nil {
				return err
			}
		}
		return nil
	}

	selector, err := sizeselect.New(sizeselect.Config{
		Model:    sizeselect.Model(*sizeModel),
		Min:      *minLen,
		Max:      *maxLen,
		ScoreMin: *scoreMin,
		ScoreMax: *scoreMax,
		Mean:     *sizeMean,
		SD:       *sizeSD,
		EdgeSD:   *sizeEdgeSD,
	})
	if err != nil {
		return err
	}

	index, err := screen.BuildCutIndexFromFASTA(*fastaPath, enzymes, digest.Options{StrictCuts: *strictCuts})
	if err != nil {
		return err
	}
	if _, err := fmt.Fprintf(stderr, "records\t%d\n", len(index.Records)); err != nil {
		return err
	}
	if _, err := fmt.Fprintf(stderr, "cached_cut_sites\t%d\n", index.CachedCutSites()); err != nil {
		return err
	}
	if _, err := fmt.Fprintf(stderr, "cache_memory_estimate_bytes\t%d\n", index.CacheMemoryEstimateBytes()); err != nil {
		return err
	}

	workers := *jobsFlag
	if workers <= 0 {
		workers = *threadsFlag
	}
	if workers <= 0 {
		workers = 1
	}
	if workers > len(pairJobs) {
		workers = len(pairJobs)
	}

	opt := digest.Options{
		AllowSame:   *allowSame,
		IncludeEnds: *includeEnds,
		StrictCuts:  *strictCuts,
	}
	return scorePairJobs(pairJobs, index, selector, opt, workers, args, stderr)
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

func buildPairJobs(enzymeNames []string, outDir string, force bool, maxPairs int) []pairJob {
	jsonDir := filepath.Join(outDir, "json")
	logDir := filepath.Join(outDir, "logs")
	jobs := make([]pairJob, 0, len(enzymeNames)*(len(enzymeNames)-1)/2)
	for i := 0; i < len(enzymeNames); i++ {
		for j := i + 1; j < len(enzymeNames); j++ {
			tag := pairTag(enzymeNames[i], enzymeNames[j])
			jsonPath := filepath.Join(jsonDir, tag+".json")
			if !force {
				if _, err := os.Stat(jsonPath); err == nil {
					continue
				}
			}
			jobs = append(jobs, pairJob{
				enzymeA: enzymeNames[i],
				enzymeB: enzymeNames[j],
				tag:     tag,
				json:    jsonPath,
				log:     filepath.Join(logDir, tag+".log"),
			})
			if maxPairs > 0 && len(jobs) >= maxPairs {
				return jobs
			}
		}
	}
	return jobs
}

func pairTag(enzymeA, enzymeB string) string {
	return safeTag.ReplaceAllString(enzymeA+"__"+enzymeB, "_")
}

func scorePairJobs(jobs []pairJob, index screen.CutIndex, selector sizeselect.Selector, opt digest.Options, workers int, command []string, stderr io.Writer) error {
	jobCh := make(chan pairJob)
	resultCh := make(chan pairResult, len(jobs))
	var wg sync.WaitGroup
	for w := 0; w < workers; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for job := range jobCh {
				resultCh <- pairResult{job: job, err: scoreOnePair(job, index, selector, opt, command)}
			}
		}()
	}
	for _, job := range jobs {
		jobCh <- job
	}
	close(jobCh)
	wg.Wait()
	close(resultCh)

	var failed []pairResult
	results := make([]pairResult, 0, len(jobs))
	for result := range resultCh {
		results = append(results, result)
		if result.err != nil {
			failed = append(failed, result)
		}
	}
	sort.Slice(results, func(i, j int) bool { return results[i].job.tag < results[j].job.tag })
	for _, result := range results {
		status := "ok"
		if result.err != nil {
			status = "error=" + result.err.Error()
		}
		if _, err := fmt.Fprintf(stderr, "%s\t%s\t%s\t%s\n", result.job.tag, status, result.job.json, result.job.log); err != nil {
			return err
		}
	}
	if len(failed) > 0 {
		return fmt.Errorf("failed_pairs\t%d", len(failed))
	}
	return nil
}

func scoreOnePair(job pairJob, index screen.CutIndex, selector sizeselect.Selector, opt digest.Options, command []string) error {
	started := time.Now()
	if err := os.MkdirAll(filepath.Dir(job.json), 0o755); err != nil {
		return err
	}
	if err := os.MkdirAll(filepath.Dir(job.log), 0o755); err != nil {
		return err
	}
	summary, err := screen.ScorePair(index, job.enzymeA, job.enzymeB, selector, opt)
	if err != nil {
		if logErr := writePairLog(job, command, started, err); logErr != nil {
			return errors.Join(err, logErr)
		}
		return err
	}
	out := pairJSON{
		SchemaVersion:   summary.SchemaVersion,
		RadigestVersion: version,
		Command:         append([]string(nil), command...),
		Enzymes:         summary.Enzymes,
		MinLength:       summary.MinLength,
		MaxLength:       summary.MaxLength,
		TotalFragments:  summary.TotalFragments,
		TotalBases:      summary.TotalBases,
		PerChromosome:   summary.PerChromosome,
		SizeSelection:   summary.SizeSelection,
		Screening:       summary.Screening,
	}
	if err := writeJSONAtomic(job.json, out); err != nil {
		if logErr := writePairLog(job, command, started, err); logErr != nil {
			return errors.Join(err, logErr)
		}
		return err
	}
	return writePairLog(job, command, started, nil)
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

func writePairLog(job pairJob, command []string, started time.Time, err error) error {
	if err := os.MkdirAll(filepath.Dir(job.log), 0o755); err != nil {
		return err
	}

	var b strings.Builder
	b.WriteString(fmt.Sprintf("pair\t%s,%s\n", job.enzymeA, job.enzymeB))
	b.WriteString(fmt.Sprintf("json\t%s\n", job.json))
	b.WriteString(fmt.Sprintf("started_at\t%s\n", started.Format(time.RFC3339Nano)))
	b.WriteString(fmt.Sprintf("finished_at\t%s\n", time.Now().Format(time.RFC3339Nano)))
	b.WriteString(fmt.Sprintf("command\t%s\n", strings.Join(command, " ")))
	if err != nil {
		b.WriteString(fmt.Sprintf("status\terror\nerror\t%s\n", err))
	} else {
		b.WriteString("status\tok\n")
	}

	return os.WriteFile(job.log, []byte(b.String()), 0o644)
}
