package digest

import (
	"sort"
	"sync"
	"strings"

	"radigest/internal/enzyme"
)

var (
	intSlicePool = sync.Pool{
		New: func() any { return make([]int, 0, 1024) },
	}
	cutSlicePool = sync.Pool{
		New: func() any { return make([]cut, 0, 1024) },
	}
)


// Fragment is half-open, 0-based [Start, End).
type Fragment struct {
	Start int
	End   int
}

type cut struct {
	pos int
	enz int
}

// Digest keeps only A-B or B-A fragments whose length ∈ [min,max].
// The first two enzymes in ens are taken as A and B.
func Digest(seq []byte, ens []enzyme.Enzyme, min, max int) []Fragment {
	if len(ens) < 2 {
		return nil
	}

	// pre-compile once
	type matcher struct {
		mask   []uint8
		offset int
	}
	m := make([]matcher, len(ens))

	for i, e := range ens {
		site := e.Recognition
		offset := e.CutIndex                              // default

		if idx := strings.IndexByte(site, '^'); idx >= 0 { // caret style
			site = site[:idx] + site[idx+1:]
			offset = idx
		} else if offset == 0 {                            // no caret + no CutIndex
			offset = len(site) / 2                        // mid-site fallback
		}

		m[i].mask = enzyme.CompileMask(site)
		m[i].offset = offset
	}

	// collect cuts per enzyme
	cutsByEnz := make([][]int, len(ens))
	for i := range cutsByEnz {                    // pooled [][]int
		cutsByEnz[i] = intSlicePool.Get().([]int)[:0]
	}
	for i, mat := range m {
		n := len(mat.mask)
		for p := 0; p <= len(seq)-n; p++ {
			if enzyme.MatchMask(mat.mask, seq[p:p+n]) {
				cutsByEnz[i] = append(cutsByEnz[i], p+mat.offset)
			}
		}
	}

	// merge → one sorted slice
	all := cutSlicePool.Get().([]cut)[:0]          // pooled []cut
	for enzID, slice := range cutsByEnz {
		for _, c := range slice {
			all = append(all, cut{pos: c, enz: enzID})
		}
	}
	sort.Slice(all, func(i, j int) bool { return all[i].pos < all[j].pos })

	// keep only AB / BA fragments
	aID, bID := 0, 1
	out := make([]Fragment, 0, len(all))
	for i := 1; i < len(all); i++ {
		l, r := all[i-1], all[i]
		if (l.enz == aID && r.enz == bID) || (l.enz == bID && r.enz == aID) {
			if ln := r.pos - l.pos; ln >= min && ln <= max {
				out = append(out, Fragment{Start: l.pos, End: r.pos})
			}
		}
	}
	// --- recycle big temps --------------------------------------------------
	for i := range cutsByEnz {
		intSlicePool.Put(cutsByEnz[i][:0])
	}
	cutSlicePool.Put(all[:0])
	
	return out
}
