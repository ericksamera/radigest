package digest

import (
	"strings"
	"sync"

	"radigest/internal/enzyme"
)

// Fragment is half-open, 0-based [Start, End).
type Fragment struct {
	Start int
	End   int
}

type matcher struct {
	mask   []uint8
	offset int
}

// Plan precompiles up to two enzymes (A,B) for fast reuse.
type Plan struct {
	m [2]matcher // m[0] = A (required), m[1] = B (optional)
}

// NewPlan compiles the first one or two enzymes.
func NewPlan(ens []enzyme.Enzyme) Plan {
	var p Plan
	n := 2
	if len(ens) < n {
		n = len(ens)
	}
	for i := 0; i < n; i++ {
		e := ens[i]
		site := e.Recognition
		offset := e.CutIndex
		if idx := strings.IndexByte(site, '^'); idx >= 0 {
			site = site[:idx] + site[idx+1:]
			offset = idx
		} else if offset == 0 {
			offset = len(site) / 2
		}
		p.m[i] = matcher{
			mask:   enzyme.CompileMask(site),
			offset: offset,
		}
	}
	return p
}

var intSlicePool = sync.Pool{
	New: func() any { return make([]int, 0, 1024) },
}

// Digest supports:
//   • single‑enzyme mode (only A configured): consecutive A cuts
//   • double‑enzyme mode (A,B): adjacent AB or BA only
func (p Plan) Digest(seq []byte, min, max int) []Fragment {
	if p.m[0].mask == nil { // no enzymes compiled
		return nil
	}

	// Scan for A and B cuts (if B is present). Slices are naturally sorted.
	aCuts := intSlicePool.Get().([]int)[:0]
	bCuts := intSlicePool.Get().([]int)[:0]
	defer func() {
		intSlicePool.Put(aCuts[:0])
		intSlicePool.Put(bCuts[:0])
	}()

	// A cuts
	{
		mat := p.m[0]
		n := len(mat.mask)
		if n > 0 && len(seq) >= n {
			for pos := 0; pos <= len(seq)-n; pos++ {
				if enzyme.MatchMask(mat.mask, seq[pos:pos+n]) {
					aCuts = append(aCuts, pos+mat.offset)
				}
			}
		}
	}

	// Single‑enzyme mode
	if p.m[1].mask == nil {
		out := make([]Fragment, 0, len(aCuts))
		for k := 1; k < len(aCuts); k++ {
			if ln := aCuts[k] - aCuts[k-1]; ln >= min && ln <= max {
				out = append(out, Fragment{Start: aCuts[k-1], End: aCuts[k]})
			}
		}
		return out
	}

	// Double‑enzyme mode: also collect B cuts
	{
		mat := p.m[1]
		n := len(mat.mask)
		if n > 0 && len(seq) >= n {
			for pos := 0; pos <= len(seq)-n; pos++ {
				if enzyme.MatchMask(mat.mask, seq[pos:pos+n]) {
					bCuts = append(bCuts, pos+mat.offset)
				}
			}
		}
	}

	// Linear two‑way merge over aCuts and bCuts → emit adjacent AB/BA only.
	out := make([]Fragment, 0, (len(aCuts)+len(bCuts))/2)
	i, j := 0, 0
	prevType := -1 // 0 = A, 1 = B
	prevPos := 0

	for i < len(aCuts) || j < len(bCuts) {
		curType, curPos := 0, 0
		if i < len(aCuts) && (j >= len(bCuts) || aCuts[i] < bCuts[j]) {
			curType, curPos = 0, aCuts[i]
			i++
		} else {
			curType, curPos = 1, bCuts[j]
			j++
		}
		if prevType != -1 && prevType != curType {
			if ln := curPos - prevPos; ln >= min && ln <= max {
				out = append(out, Fragment{Start: prevPos, End: curPos})
			}
		}
		prevType, prevPos = curType, curPos
	}
	return out
}

// Back‑compat convenience: compile plan per call.
func Digest(seq []byte, ens []enzyme.Enzyme, min, max int) []Fragment {
	return NewPlan(ens).Digest(seq, min, max)
}
