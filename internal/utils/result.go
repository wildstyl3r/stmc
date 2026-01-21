package utils

import "reflect"

type CoreResult struct {
	ModelName                      string  `csv:"Model name"`
	CathodeFallLength              float64 `csv:"Sheath length" units:"Length:1"`
	PressureCathodeFallLength      float64 `csv:"pd" units:"Length:1"`
	ReducedFieldAtCathode          float64 `csv:"E/n at cathode"`
	ReducedFieldAtSheathCenter     float64 `csv:"E/n at mid-sheath"`
	GlobalMeanFreePath             float64 `csv:"Global mean free path" units:"Length:1"`
	Voltage                        float64 `csv:"Voltage" units:"Voltage:1"`
	Pressure                       float64 `csv:"p" units:"Pressure:1"`
	ElectronsReturned              float64 `csv:"electrons returned"`
	ElectronsReturnedMargin        float64 `csv:"electrons returned margin"`
	IonizingCathodeElectrons       float64 `csv:"ionizing cathode electrons"`
	IonizingCathodeElectronsMargin float64 `csv:"ionizing cathode electrons margin"`
}

type ResultInterface interface {
	Index() float64
	SetModelName(string)
}

func (row *CoreResult) Index() float64 {
	return row.ReducedFieldAtCathode
}

func (row *CoreResult) SetModelName(s string) {
	row.ModelName = s
}

func HomogeneousResultInterfaceSliceToStructSlice(s []ResultInterface) any {
	if len(s) == 0 {
		return nil
	}
	t := reflect.TypeOf(s[0])
	sliceType := reflect.SliceOf(t)
	structSlice := reflect.MakeSlice(sliceType, len(s), len(s))
	for i := range s {
		structSlice.Index(i).Set(reflect.ValueOf(s[i]))
	}
	return structSlice.Interface()
}
