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
