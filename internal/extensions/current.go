package extensions

import (
	"math"

	"github.com/wildstyl3r/stmc/internal/config"
	"github.com/wildstyl3r/stmc/internal/constants"
	"github.com/wildstyl3r/stmc/internal/model"
	"github.com/wildstyl3r/stmc/internal/utils"
)

func CurrentDensityEquation(parameters config.ModelParameters) float64 {
	dc, secondaryEmissionCoefficient, Vc, N := parameters.CathodeFallLength, gammaAnalyticF(parameters), parameters.CathodeFallPotential, parameters.GasDensity
	ionDriftVelocity := utils.IonDriftVelocity[parameters.Species]
	return (1 + secondaryEmissionCoefficient) * ionDriftVelocity(Vc, dc, N) * 2 * Vc / (dc * dc) * constants.FreeSpacePermittivityE0
}

func CurrentDensityRelativeError(sourceIntegral, sourceIntegralMargin, dc, dcMargin float64) float64 {
	return sourceIntegralMargin/(sourceIntegral*(sourceIntegral+1)) + 2.5*dcMargin/dc
}

func CurrentDensityCalculation(parameters config.ModelParameters, outputDir, modelName string) (dataRow CurrentDensityDataRow, altRow *CurrentDensityDataRow, finalModel, altModel *model.Model) {
	minDc := parameters.CathodeFallLengthPrecision
	maxDc := parameters.SimulationLength()
	numberOfSteps := int((maxDc - minDc) / parameters.CathodeFallLengthPrecision)
	return GeneralizedCalculation(parameters, numberOfSteps, minDc, parameters.CathodeFallLengthPrecision,
		func(dc float64) config.ModelParameters {
			newParameters := parameters
			newParameters.CathodeFallLength = dc
			return newParameters
		},
		func(optResult OptimizationResult, variance float64, m *model.Model) CurrentDensityDataRow {
			optResult.CathodeCurrentDensity = CurrentDensityEquation(m.Parameters)
			optResult.CathodeCurrent = optResult.CathodeCurrentDensity * (math.Pi * parameters.CathodeRadius * parameters.CathodeRadius)
			optResult.CathodeCurrentDensityPerPressureSquared = optResult.CathodeCurrentDensity / (optResult.Pressure * optResult.Pressure)
			return CurrentDensityDataRow{
				OptimizationResult:                            optResult,
				CathodeCurrentDensityMargin:                   optResult.CathodeCurrentDensity * CurrentDensityRelativeError(optResult.SourceIntegralMonteCarlo, math.Abs(optResult.SourceIntegralDifference)+optResult.SourceIntegralMargin, optResult.CathodeFallLength, optResult.CathodeFallLengthMargin),
				CathodeCurrentDensityPerPressureSquaredMargin: optResult.CathodeCurrentDensity / (parameters.Pressure * parameters.Pressure) * CurrentDensityRelativeError(optResult.SourceIntegralMonteCarlo, math.Abs(optResult.SourceIntegralDifference)+optResult.SourceIntegralMargin, optResult.CathodeFallLength, optResult.CathodeFallLengthMargin),
			}
		}, func(mp config.ModelParameters, simc, sia float64) bool {
			return (sia > 0 && simc > 2*sia && simc > 30.)
		}, true, outputDir, modelName)
}

func CurrentDensityCalculation2(parameters config.ModelParameters, outputDir, modelName string) (dataRow CurrentDensityDataRow, altRow *CurrentDensityDataRow, finalModel, altModel *model.Model) {
	minDc := parameters.CathodeFallLengthPrecision
	maxDc := parameters.SimulationLength()
	numberOfSteps := int((maxDc - minDc) / parameters.CathodeFallLengthPrecision)
	return GeneralizedCalculation(parameters, numberOfSteps, minDc, parameters.CathodeFallLengthPrecision,
		func(dc float64) config.ModelParameters {
			newParameters := parameters
			newParameters.CathodeFallLength = dc
			return newParameters
		},
		func(optResult OptimizationResult, variance float64, m *model.Model) CurrentDensityDataRow {
			optResult.CathodeCurrentDensity = CurrentDensityEquation(m.Parameters)
			optResult.CathodeCurrentDensityPerPressureSquared = optResult.CathodeCurrentDensity / (optResult.Pressure * optResult.Pressure)
			currentDensityMargin := optResult.CathodeCurrentDensity * CurrentDensityRelativeError(optResult.SourceIntegralMonteCarlo, math.Abs(optResult.SourceIntegralDifference)+optResult.SourceIntegralMargin, optResult.CathodeFallLength, optResult.CathodeFallLengthMargin)
			return CurrentDensityDataRow{
				OptimizationResult:                            optResult,
				CathodeCurrentMargin:                          currentDensityMargin * (math.Pi * parameters.CathodeRadius * parameters.CathodeRadius),
				CathodeCurrentDensityMargin:                   currentDensityMargin,
				CathodeCurrentDensityPerPressureSquaredMargin: currentDensityMargin / (parameters.Pressure * parameters.Pressure),
			}
		}, func(mp config.ModelParameters, simc, sia float64) bool {
			return (sia > 0 && simc > 2*sia && simc > 30.)
		}, true, outputDir, modelName)
}

type CurrentDensityDataRow struct {
	OptimizationResult
	CathodeCurrentMargin                          float64 `csv:"I margin" units:"Current:1"`
	CathodeCurrentDensityMargin                   float64 `csv:"j margin" units:"Current:1,Length:-2"`
	CathodeCurrentDensityPerPressureSquaredMargin float64 `csv:"j/p2 margin" units:"Current:1,Length:-2,Pressure:-2"`
}
