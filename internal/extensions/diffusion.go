package extensions

import (
	"math"

	"github.com/wildstyl3r/stmc/internal/model"
)

const (
	CharacteristicDiffusionScaleKey  = "CharacteristicDiffusionScale"
	AmbipolarDiffusionCoefficientKey = "AmbipolarDiffusionCoefficient"
)

func DiffusionConstants(m *model.Model) ([]string, []any, error) {
	var cylindricalCharacteristicDiffusionLength float64
	if m.Parameters.SimplifiedDiffusionScale {
		cylindricalCharacteristicDiffusionLength = m.Parameters.TubeRadius / 2.405
	} else {
		var L float64
		if m.Parameters.ParallelPlaneHollowCathode {
			L = m.Parameters.GapLength * 2 // was halved at initialization
		} else {
			L = m.Parameters.GapLength
		}
		cylindricalCharacteristicDiffusionLength = L * m.Parameters.TubeRadius /
			math.Sqrt(math.Pi*math.Pi*m.Parameters.TubeRadius*m.Parameters.TubeRadius+2.405*2.405*L*L)
	}
	ambipolarDiffusionCoefficient := m.Parameters.IonMobility * m.Parameters.SlowElectronTemperature
	return []string{CharacteristicDiffusionScaleKey, AmbipolarDiffusionCoefficientKey}, []any{cylindricalCharacteristicDiffusionLength, ambipolarDiffusionCoefficient}, nil
}
