package extensions

import (
	"github.com/wildstyl3r/stmc/internal/config"
	"github.com/wildstyl3r/stmc/internal/constants"
	"github.com/wildstyl3r/stmc/internal/utils"
)

// not applicable for UniformField
func gammaAnalyticPhelpsF(parameters *config.ModelParameters) float64 {
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
		g := gammaAnalyticF(&parameters)
		return 0. < g
	}, lowerBound, upperBound, parameters.CathodeFallLengthPrecision*0.01)

	_, to = utils.BinarySearch(func(dc float64) bool {
		parameters.CathodeFallLength = dc
		g := gammaAnalyticF(&parameters)
		return 1 < g // might be false everywhere in the gap, but as close as possible to true domain
	}, lowerBound, upperBound, parameters.CathodeFallLengthPrecision*0.01)
	return from, to
}

type OptimizationResult struct {
	utils.SheathResult
	CathodeFallLengthMargin                 float64 `csv:"Sheath length margin" units:"Length:1"`
	PressureCathodeFallLengthMargin         float64 `csv:"pd margin" units:"Length:1"`
	SourceIntegralDifference                float64 `csv:"Source integral loss"`
	SourceIntegralAnalytic                  float64 `csv:"Source integral analytic"`
	SourceIntegralMonteCarlo                float64 `csv:"Source integral Monte Carlo"`
	SourceIntegralMargin                    float64 `csv:"Source integral margin"`
	EffectiveGammaAnalytic                  float64 `csv:"effective $\\gamma$ analytic"`
	EffectiveGammaMonteCarlo                float64 `csv:"effective $\\gamma$ Monte Carlo"`
	SurfaceGammaAnalytic                    float64 `csv:"surface $\\gamma$ analytic"`
	SurfaceGammaMonteCarlo                  float64 `csv:"surface $\\gamma$ Monte Carlo"`
	MeanElectronEnergyAtAnode               float64 `csv:"Mean electron energy at anode" units:"Energy:1"`
	MeanFreePathAnode                       float64 `csv:"Mean free path at anode" units:"Length:1"`
	CathodeCurrent                          float64 `csv:"I" units:"Current:1"`
	CathodeCurrentDensity                   float64 `csv:"j" units:"Current:1,Length:-2"`
	CathodeCurrentDensityPerPressureSquared float64 `csv:"j/p2" units:"Current:1,Length:-2,Pressure:-2"`
}
