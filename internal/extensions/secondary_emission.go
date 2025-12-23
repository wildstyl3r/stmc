package extensions

import (
	"github.com/wildstyl3r/stmc/internal/config"
	"github.com/wildstyl3r/stmc/internal/constants"
	"github.com/wildstyl3r/stmc/internal/utils"
)

// not applicable for UniformField
func gammaAnalyticPhelpsF(parameters config.ModelParameters) float64 {
	dc, j, Vc, N := parameters.CathodeFallLength, parameters.CathodeCurrentDensity, parameters.CathodeFallPotential, parameters.GasDensity
	ionDriftVelocity := utils.IonDriftVelocity[parameters.GetSpeciesString()]
	return j*dc*dc/(2.*Vc*ionDriftVelocity(Vc, dc, N)*constants.FreeSpacePermittivityE0) - 1. //-> 0
}

func EstimateCathodeFallLengthLimits(parameters config.ModelParameters) (from float64, to float64) {
	upperBound := parameters.SimulationLength()
	upperBound -= parameters.CathodeFallLengthPrecision * 0.01

	lowerBound := 1e-9

	from, _ = utils.BinarySearch(func(dc float64) bool {
		parameters.CathodeFallLength = dc
		g := gammaAnalyticF(parameters)
		return 0. < g
	}, lowerBound, upperBound, parameters.CathodeFallLengthPrecision*0.01)

	_, to = utils.BinarySearch(func(dc float64) bool {
		parameters.CathodeFallLength = dc
		g := gammaAnalyticF(parameters)
		return 1 < g // might be false everywhere in the gap, but as close as possible to true domain
	}, lowerBound, upperBound, parameters.CathodeFallLengthPrecision*0.01)
	return from, to
}

type OptimizationResult struct {
	ModelName                               string  `csv:"Model name"`
	ReducedFieldAtCathode                   float64 `csv:"E/n at cathode"`
	ReducedFieldAtSheathCenter              float64 `csv:"E/n at mid-sheath"`
	CathodeFallLength                       float64 `csv:"Sheath length" units:"Length:1"`
	CathodeFallLengthMargin                 float64 `csv:"Sheath length margin" units:"Length:1"`
	PressureCathodeFallLength               float64 `csv:"pd" units:"Length:1"`
	PressureCathodeFallLengthMargin         float64 `csv:"pd margin" units:"Length:1"`
	SourceIntegralDifference                float64 `csv:"Source integral loss"`
	SourceIntegralAnalytic                  float64 `csv:"Source integral analytic"`
	SourceIntegralMonteCarlo                float64 `csv:"Source integral Monte Carlo"`
	SourceIntegralMargin                    float64 `csv:"Source integral margin"`
	GammaAnalytic                           float64 `csv:"$\\gamma$ analytic"`
	GammaMonteCarlo                         float64 `csv:"$\\gamma$ Monte Carlo"`
	MeanElectronEnergyAtAnode               float64 `csv:"Mean electron energy at anode" units:"Energy:1"`
	MeanFreePathAnode                       float64 `csv:"Mean free path at anode" units:"Length:1"`
	GlobalMeanFreePath                      float64 `csv:"Global mean free path" units:"Length:1"`
	Voltage                                 float64 `csv:"Voltage" units:"Voltage:1"`
	Pressure                                float64 `csv:"p" units:"Pressure:1"`
	CathodeCurrent                          float64 `csv:"I" units:"Current:1"`
	CathodeCurrentDensity                   float64 `csv:"j" units:"Current:1,Length:-2"`
	CathodeCurrentDensityPerPressureSquared float64 `csv:"j/p2" units:"Current:1,Length:-2,Pressure:-2"`
}

func (row OptimizationResult) Index() float64 {
	return row.ReducedFieldAtCathode
}
