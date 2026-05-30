package digest

import (
	"bytes"
	"fmt"
	"strings"

	"github.com/ericksamera/radigest/internal/enzyme"
)

// Fragment is half-open, 0-based [Start, End).
type Fragment struct {
	Start int
	End   int
}

// Stats summarizes kept digest fragments without materializing Fragment values.
type Stats struct {
	Fragments int
	Bases     int
}

func (s *Stats) addIfKept(start, end, min, max int) {
	if ln := end - start; ln >= min && ln <= max {
		s.Fragments++
		s.Bases += ln
	}
}

func (s *Stats) addTerminalIfKept(start, end, min, max int) {
	if end <= start {
		return
	}
	s.addIfKept(start, end, min, max)
}

type matcher struct {
	exact  []byte
	mask   []uint8
	anchor int
	offset int
}

type Options struct {
	AllowSame   bool // keep AA/BB neighbors in double digest
	StrictCuts  bool // error if site has no caret and CutIndex==0 (mid-site fallback)
	IncludeEnds bool // also emit terminal chromosome/contig-end fragments
}

// Plan precompiles up to two enzymes (A,B) for fast reuse.
type Plan struct {
	m           [2]matcher // m[0] = A (required), m[1] = B (optional)
	allowSame   bool
	includeEnds bool
}

func NewPlanWithOptions(ens []enzyme.Enzyme, opt Options) Plan {
	p, err := TryNewPlanWithOptions(ens, opt)
	if err != nil {
		panic(err)
	}
	return p
}

func TryNewPlanWithOptions(ens []enzyme.Enzyme, opt Options) (Plan, error) {
	var p Plan
	p.allowSame = opt.AllowSame
	p.includeEnds = opt.IncludeEnds

	n := 2
	if len(ens) < n {
		n = len(ens)
	}
	for i := 0; i < n; i++ {
		e := ens[i]
		site := e.Recognition
		offset := e.CutIndex
		usedFallback := false
		if idx := strings.IndexByte(site, '^'); idx >= 0 {
			site = site[:idx] + site[idx+1:]
			offset = idx
		} else if offset == 0 {
			usedFallback = true
			offset = len(site) / 2
		}
		if site == "" {
			return Plan{}, fmt.Errorf("enzyme %s: empty recognition site", e.Name)
		}
		if offset < 0 || offset > len(site) {
			return Plan{}, fmt.Errorf("enzyme %s: cut offset %d outside recognition site length %d", e.Name, offset, len(site))
		}
		if opt.StrictCuts && usedFallback {
			return Plan{}, fmt.Errorf("enzyme %s: no caret and CutIndex==0 (mid-site fallback disabled by -strict-cuts)", e.Name)
		}
		mask, err := enzyme.CompileMaskChecked(site)
		if err != nil {
			return Plan{}, fmt.Errorf("enzyme %s recognition %q: %w", e.Name, e.Recognition, err)
		}
		mat := matcher{
			mask:   mask,
			anchor: enzyme.BestMaskAnchor(mask),
			offset: offset,
		}
		if enzyme.IsExactACGT(site) {
			mat.exact = []byte(strings.ToUpper(site))
		}
		p.m[i] = mat
	}
	return p, nil
}

// Back-compat.
func NewPlan(ens []enzyme.Enzyme) Plan { return NewPlanWithOptions(ens, Options{}) }

type cutScanner struct {
	mat matcher
	seq []byte
	pos int
}

func newCutScanner(mat matcher, seq []byte) cutScanner {
	return cutScanner{mat: mat, seq: seq}
}

func (s *cutScanner) next() (int, bool) {
	if len(s.mat.exact) > 0 {
		return s.nextExact()
	}
	return s.nextMask()
}

func (s *cutScanner) nextMask() (int, bool) {
	n := len(s.mat.mask)
	if n == 0 || len(s.seq) < n {
		return 0, false
	}
	for s.pos <= len(s.seq)-n {
		pos := s.pos
		s.pos++
		if enzyme.MatchMaskAt(s.mat.mask, s.mat.anchor, s.seq[pos:pos+n]) {
			return pos + s.mat.offset, true
		}
	}
	return 0, false
}

func (s *cutScanner) nextExact() (int, bool) {
	n := len(s.mat.exact)
	if n == 0 || len(s.seq) < n || s.pos > len(s.seq)-n {
		return 0, false
	}

	idx := bytes.Index(s.seq[s.pos:], s.mat.exact)
	if idx < 0 {
		return 0, false
	}

	siteStart := s.pos + idx
	s.pos = siteStart + 1 // preserve overlapping motif detection
	return siteStart + s.mat.offset, true
}

func emitIfKept(start, end, min, max int, emit func(Fragment) error) error {
	if ln := end - start; ln >= min && ln <= max {
		return emit(Fragment{Start: start, End: end})
	}
	return nil
}

func emitTerminalIfKept(start, end, min, max int, emit func(Fragment) error) error {
	if end <= start {
		return nil
	}
	return emitIfKept(start, end, min, max, emit)
}

// CutsEach streams sorted cut coordinates for the first enzyme in the plan.
// Cut coordinates are motif start plus cut offset. The callback is invoked in
// deterministic genomic cut-coordinate order. If emit returns an error,
// scanning stops and that error is returned.
func (p Plan) CutsEach(seq []byte, emit func(int) error) error {
	if p.m[0].mask == nil { // no enzymes compiled
		return nil
	}
	if emit == nil {
		return fmt.Errorf("digest cut emit callback is nil")
	}

	scan := newCutScanner(p.m[0], seq)
	for {
		cut, ok := scan.next()
		if !ok {
			return nil
		}
		if err := emit(cut); err != nil {
			return err
		}
	}
}

// Cuts returns sorted cut coordinates for the first enzyme in the plan.
func (p Plan) Cuts(seq []byte) []int {
	if p.m[0].mask == nil {
		return nil
	}
	cuts := make([]int, 0)
	_ = p.CutsEach(seq, func(cut int) error {
		cuts = append(cuts, cut)
		return nil
	})
	return cuts
}

// DigestEach streams kept fragments to emit without materializing cut arrays or
// a per-chromosome []Fragment. It supports the same modes as Digest:
//   - single-enzyme mode (only A configured): consecutive A cuts
//   - double-enzyme mode (A,B): adjacent AB/BA only, or AA/BB too if AllowSame
//   - optional terminal chromosome/contig-end fragments if IncludeEnds is set
//
// The callback is invoked in deterministic genomic cut-coordinate order. If emit
// returns an error, scanning stops and that error is returned.
func (p Plan) DigestEach(seq []byte, min, max int, emit func(Fragment) error) error {
	if p.m[0].mask == nil { // no enzymes compiled
		return nil
	}
	if emit == nil {
		return fmt.Errorf("digest emit callback is nil")
	}

	aScan := newCutScanner(p.m[0], seq)
	aPos, aOK := aScan.next()

	// Single-enzyme mode: only the previous cut coordinate is needed.
	if p.m[1].mask == nil {
		if !aOK {
			if p.includeEnds {
				return emitTerminalIfKept(0, len(seq), min, max, emit)
			}
			return nil
		}
		if p.includeEnds {
			if err := emitTerminalIfKept(0, aPos, min, max, emit); err != nil {
				return err
			}
		}
		prevPos := aPos
		for {
			pos, ok := aScan.next()
			if !ok {
				if p.includeEnds {
					return emitTerminalIfKept(prevPos, len(seq), min, max, emit)
				}
				return nil
			}
			if err := emitIfKept(prevPos, pos, min, max, emit); err != nil {
				return err
			}
			prevPos = pos
		}
	}

	// Double-enzyme mode: merge the two naturally sorted cut-coordinate streams.
	bScan := newCutScanner(p.m[1], seq)
	bPos, bOK := bScan.next()
	prevType := -1 // 0=A, 1=B
	prevPos := 0
	sawCut := false
	lastPos := 0

	for aOK || bOK {
		var pos int
		if aOK && (!bOK || aPos <= bPos) {
			pos = aPos
		} else {
			pos = bPos
		}
		hasA := aOK && aPos == pos
		hasB := bOK && bPos == pos

		if p.includeEnds && !sawCut {
			if err := emitTerminalIfKept(0, pos, min, max, emit); err != nil {
				return err
			}
		}
		sawCut = true
		lastPos = pos

		if hasA && hasB {
			// Coincident cuts are barriers. Emit one zero-length fragment for the
			// site if the caller's size range allows it, then reset adjacency so no
			// fragment bridges across the coincident cut.
			if err := emitIfKept(pos, pos, min, max, emit); err != nil {
				return err
			}
			aPos, aOK = aScan.next()
			bPos, bOK = bScan.next()
			prevType = -1
			prevPos = pos
			continue
		}

		curType := 0
		if hasA {
			aPos, aOK = aScan.next()
		} else {
			curType = 1
			bPos, bOK = bScan.next()
		}

		if prevType != -1 && (p.allowSame || prevType != curType) {
			if err := emitIfKept(prevPos, pos, min, max, emit); err != nil {
				return err
			}
		}
		prevType, prevPos = curType, pos
	}
	if p.includeEnds {
		if !sawCut {
			return emitTerminalIfKept(0, len(seq), min, max, emit)
		}
		return emitTerminalIfKept(lastPos, len(seq), min, max, emit)
	}
	return nil
}

// DigestCutsEach streams kept fragments from precomputed sorted cut-coordinate
// slices. If cutsB is nil, single-enzyme semantics are used for cutsA. If cutsB
// is non-nil, double-digest semantics are used for cutsA and cutsB.
//
// The behavior matches Plan.DigestEach for single digest, double digest,
// coincident cuts, AA/BB suppression, AllowSame, and IncludeEnds, assuming the
// input cut slices are sorted in ascending cut-coordinate order.
func DigestCutsEach(cutsA, cutsB []int, seqLen, min, max int, opt Options, emit func(Fragment) error) error {
	if emit == nil {
		return fmt.Errorf("digest emit callback is nil")
	}
	if seqLen < 0 {
		return fmt.Errorf("digest sequence length is negative: %d", seqLen)
	}
	if cutsB == nil {
		return digestSingleCutsEach(cutsA, seqLen, min, max, opt, emit)
	}
	return digestDoubleCutsEach(cutsA, cutsB, seqLen, min, max, opt, emit)
}

func digestSingleCutsEach(cuts []int, seqLen, min, max int, opt Options, emit func(Fragment) error) error {
	if len(cuts) == 0 {
		if opt.IncludeEnds {
			return emitTerminalIfKept(0, seqLen, min, max, emit)
		}
		return nil
	}
	if opt.IncludeEnds {
		if err := emitTerminalIfKept(0, cuts[0], min, max, emit); err != nil {
			return err
		}
	}
	prevPos := cuts[0]
	for _, pos := range cuts[1:] {
		if err := emitIfKept(prevPos, pos, min, max, emit); err != nil {
			return err
		}
		prevPos = pos
	}
	if opt.IncludeEnds {
		return emitTerminalIfKept(prevPos, seqLen, min, max, emit)
	}
	return nil
}

func digestDoubleCutsEach(cutsA, cutsB []int, seqLen, min, max int, opt Options, emit func(Fragment) error) error {
	ai, bi := 0, 0
	prevType := -1 // 0=A, 1=B
	prevPos := 0
	sawCut := false
	lastPos := 0

	for ai < len(cutsA) || bi < len(cutsB) {
		var pos int
		if ai < len(cutsA) && (bi >= len(cutsB) || cutsA[ai] <= cutsB[bi]) {
			pos = cutsA[ai]
		} else {
			pos = cutsB[bi]
		}
		hasA := ai < len(cutsA) && cutsA[ai] == pos
		hasB := bi < len(cutsB) && cutsB[bi] == pos

		if opt.IncludeEnds && !sawCut {
			if err := emitTerminalIfKept(0, pos, min, max, emit); err != nil {
				return err
			}
		}
		sawCut = true
		lastPos = pos

		if hasA && hasB {
			// Coincident cuts are barriers. Emit one zero-length fragment for the
			// site if the caller's size range allows it, then reset adjacency so no
			// fragment bridges across the coincident cut.
			if err := emitIfKept(pos, pos, min, max, emit); err != nil {
				return err
			}
			ai++
			bi++
			prevType = -1
			prevPos = pos
			continue
		}

		curType := 0
		if hasA {
			ai++
		} else {
			curType = 1
			bi++
		}

		if prevType != -1 && (opt.AllowSame || prevType != curType) {
			if err := emitIfKept(prevPos, pos, min, max, emit); err != nil {
				return err
			}
		}
		prevType, prevPos = curType, pos
	}
	if opt.IncludeEnds {
		if !sawCut {
			return emitTerminalIfKept(0, seqLen, min, max, emit)
		}
		return emitTerminalIfKept(lastPos, seqLen, min, max, emit)
	}
	return nil
}

// DigestCuts returns kept fragments from precomputed sorted cut-coordinate
// slices. A nil cutsB selects single-enzyme semantics; a non-nil cutsB selects
// double-digest semantics even when the second enzyme has zero cuts.
func DigestCuts(cutsA, cutsB []int, seqLen, min, max int, opt Options) []Fragment {
	out := make([]Fragment, 0)
	_ = DigestCutsEach(cutsA, cutsB, seqLen, min, max, opt, func(fr Fragment) error {
		out = append(out, fr)
		return nil
	})
	return out
}

// DigestStats returns hard-window fragment counts and bases without constructing
// Fragment values or invoking a per-fragment callback. It uses the same fragment
// semantics as DigestEach for single/double digest, coincident cuts, same-enzyme
// adjacency, and optional terminal fragments.
func (p Plan) DigestStats(seq []byte, min, max int) Stats {
	var stats Stats
	if p.m[0].mask == nil { // no enzymes compiled
		return stats
	}

	aScan := newCutScanner(p.m[0], seq)
	aPos, aOK := aScan.next()

	// Single-enzyme mode: only the previous cut coordinate is needed.
	if p.m[1].mask == nil {
		if !aOK {
			if p.includeEnds {
				stats.addTerminalIfKept(0, len(seq), min, max)
			}
			return stats
		}
		if p.includeEnds {
			stats.addTerminalIfKept(0, aPos, min, max)
		}
		prevPos := aPos
		for {
			pos, ok := aScan.next()
			if !ok {
				if p.includeEnds {
					stats.addTerminalIfKept(prevPos, len(seq), min, max)
				}
				return stats
			}
			stats.addIfKept(prevPos, pos, min, max)
			prevPos = pos
		}
	}

	// Double-enzyme mode: merge the two naturally sorted cut-coordinate streams.
	bScan := newCutScanner(p.m[1], seq)
	bPos, bOK := bScan.next()
	prevType := -1 // 0=A, 1=B
	prevPos := 0
	sawCut := false
	lastPos := 0

	for aOK || bOK {
		var pos int
		if aOK && (!bOK || aPos <= bPos) {
			pos = aPos
		} else {
			pos = bPos
		}
		hasA := aOK && aPos == pos
		hasB := bOK && bPos == pos

		if p.includeEnds && !sawCut {
			stats.addTerminalIfKept(0, pos, min, max)
		}
		sawCut = true
		lastPos = pos

		if hasA && hasB {
			// Coincident cuts are barriers. Count one zero-length fragment for
			// the site if the caller's size range allows it, then reset
			// adjacency so no fragment bridges across the coincident cut.
			stats.addIfKept(pos, pos, min, max)
			aPos, aOK = aScan.next()
			bPos, bOK = bScan.next()
			prevType = -1
			prevPos = pos
			continue
		}

		curType := 0
		if hasA {
			aPos, aOK = aScan.next()
		} else {
			curType = 1
			bPos, bOK = bScan.next()
		}

		if prevType != -1 && (p.allowSame || prevType != curType) {
			stats.addIfKept(prevPos, pos, min, max)
		}
		prevType, prevPos = curType, pos
	}
	if p.includeEnds {
		if !sawCut {
			stats.addTerminalIfKept(0, len(seq), min, max)
			return stats
		}
		stats.addTerminalIfKept(lastPos, len(seq), min, max)
	}
	return stats
}

// Digest supports:
//   - single-enzyme mode (only A configured): consecutive A cuts
//   - double-enzyme mode (A,B): adjacent AB/BA only (or AA/BB too if allowSame)
func (p Plan) Digest(seq []byte, min, max int) []Fragment {
	if p.m[0].mask == nil { // no enzymes compiled
		return nil
	}
	out := make([]Fragment, 0)
	_ = p.DigestEach(seq, min, max, func(fr Fragment) error {
		out = append(out, fr)
		return nil
	})
	return out
}

// Back-compat convenience: compile plan per call.
func Digest(seq []byte, ens []enzyme.Enzyme, min, max int) []Fragment {
	return NewPlan(ens).Digest(seq, min, max)
}
