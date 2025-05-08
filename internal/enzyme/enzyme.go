// internal/enzyme/enzyme.go
//go:generate go run ./internal/enzyme/cmd/gen -in internal/enzyme/enzymes.json -out internal/enzyme/enzymes_generated.go
package enzyme

type Enzyme struct {
    Name        string
    Recognition string
    CutIndex    int // 0‑based offset from start of site
}

// dummy DB – will be generated later
// var DB = map[string]Enzyme{
//     "EcoRI": {Name: "EcoRI", Recognition: "G^AATTC"},
//     "MseI":  {Name: "MseI", Recognition: "T^TAA"},
// }

func Get(name string) (Enzyme, bool) {
    e, ok := DB[name]
    return e, ok
}
