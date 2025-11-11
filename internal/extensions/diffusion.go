package extensions

import (
	"math"

	"github.com/wildstyl3r/stmc/internal/constants"
	"github.com/wildstyl3r/stmc/internal/model"
	"github.com/wildstyl3r/stmc/internal/utils"
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
		L := m.Parameters.GapLength
		if math.IsInf(m.Parameters.TubeRadius, 1) {
			cylindricalCharacteristicDiffusionLength = L / math.Pi
		} else {
			cylindricalCharacteristicDiffusionLength = L * m.Parameters.TubeRadius /
				math.Sqrt(math.Pi*math.Pi*m.Parameters.TubeRadius*m.Parameters.TubeRadius+2.405*2.405*L*L)
		}

	}
	ambipolarDiffusionCoefficient := 2. / 3. * utils.WeakFieldConstantIonMobilityPerTorr[m.Parameters.Species] * (m.Parameters.Pressure * constants.Torr) * m.Parameters.SlowElectronTemperature
	return []string{CharacteristicDiffusionScaleKey, AmbipolarDiffusionCoefficientKey}, []any{cylindricalCharacteristicDiffusionLength, ambipolarDiffusionCoefficient}, nil
}
