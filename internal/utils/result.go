package utils

import "cmp"

type CoreResult struct {
	ModelName                  string  `csv:"Model name"`
	CathodeFallLength          float64 `csv:"Sheath length" units:"Length:1"`
	PressureCathodeFallLength  float64 `csv:"pd" units:"Length:1"`
	ReducedFieldAtCathode      float64 `csv:"E/n at cathode"`
	ReducedFieldAtSheathCenter float64 `csv:"E/n at mid-sheath"`
	GlobalMeanFreePath         float64 `csv:"Global mean free path" units:"Length:1"`
	Voltage                    float64 `csv:"Voltage" units:"Voltage:1"`
	Pressure                   float64 `csv:"p" units:"Pressure:1"`
	ElectronsReturned          float64 `csv:"electrons returned"`
	ElectronsReturnedMargin    float64 `csv:"electrons returned margin"`
}

type Indexable[T cmp.Ordered] interface {
	Index() T
}

func (row CoreResult) Index() float64 {
	return row.ReducedFieldAtCathode
}
