package utils

type SheathResult struct {
	Result
	CathodeFallLength          float64 `csv:"Sheath length" units:"Length:1"`
	PressureCathodeFallLength  float64 `csv:"pd" units:"Length:1"`
	ReducedFieldAtSheathCenter float64 `csv:"E/n at mid-sheath"`
}
