package sim

import (
	"math/rand"
	"time"
)

// Make returns an upper‑case DNA sequence of given length with ~gc fraction GC.
// If seed==0 we use a time-based seed; otherwise results are reproducible.
func Make(length int, gc float64, seed int64) []byte {
	if length <= 0 {
		return []byte{}
	}
	if gc < 0 {
		gc = 0
	}
	if gc > 1 {
		gc = 1
	}
	if seed == 0 {
		seed = time.Now().UnixNano()
	}
	r := rand.New(rand.NewSource(seed))

	gcCount := int(float64(length)*gc + 0.5) // nearest integer
	if gcCount < 0 {
		gcCount = 0
	}
	if gcCount > length {
		gcCount = length
	}
	atCount := length - gcCount

	seq := make([]byte, length)

	// Fill exact composition.
	for i := 0; i < gcCount; i++ {
		if r.Intn(2) == 0 {
			seq[i] = 'G'
		} else {
			seq[i] = 'C'
		}
	}
	for i := gcCount; i < gcCount+atCount; i++ {
		if r.Intn(2) == 0 {
			seq[i] = 'A'
		} else {
			seq[i] = 'T'
		}
	}

	// Shuffle to disperse bases.
	for i := length - 1; i > 0; i-- {
		j := r.Intn(i + 1)
		seq[i], seq[j] = seq[j], seq[i]
	}
	return seq
}
