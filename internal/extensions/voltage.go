package extensions

import (
	"math"

	"github.com/wildstyl3r/stmc/internal/config"
	"github.com/wildstyl3r/stmc/internal/constants"
	"github.com/wildstyl3r/stmc/internal/messages"
	"github.com/wildstyl3r/stmc/internal/model"
	"github.com/wildstyl3r/stmc/internal/utils"
)

func VoltageFromDc(parameters config.ModelParameters, dc float64) float64 {
	parameters.CathodeFallLength = dc
	_, maxV := utils.BinarySearch(func(V float64) bool {
		parameters.CathodeFallPotential = V
		phelps := sourceIntegralAnalyticPhelpsF(&parameters)
		return phelps > 0
	}, 1500, 10, 0.000001)
	lV, rV := utils.BinarySearch(func(V float64) bool {
		parameters.CathodeFallPotential = V
		phelps := sourceIntegralAnalyticPhelpsF(&parameters)
		other := sourceIntegralAnalyticF(&parameters)
		return phelps-other > 0.
	}, 1, maxV, 0.000001)
	return 0.5 * (lV + rV)
}

func MaxVoltage(parameters config.ModelParameters) float64 {
	return VoltageFromDc(parameters, parameters.SimulationLength()-parameters.CathodeFallLengthPrecision)
}

func VoltageFromDcSimplified(currentDensity, secondaryEmission, dc, gasDensity float64, species string) float64 {
	return math.Pow(currentDensity*math.Pow(dc, 2.5)*math.Sqrt(constants.Townsend*gasDensity)/(2*(1+secondaryEmission)*utils.SimplifiedIonDriftVelocityCoefficient[species]*constants.FreeSpacePermittivityE0), 2./3.)
}

func DcFromVoltage(parameters config.ModelParameters, voltage float64) float64 {
	parameters.CathodeFallPotential = voltage
	_, minDc := utils.BinarySearch(func(dc float64) bool {
		parameters.CathodeFallLength = dc
		return sourceIntegralAnalyticPhelpsF(&parameters) > 0
	}, parameters.CathodeFallLengthPrecision, parameters.SimulationLength(), parameters.CathodeFallLengthPrecision*0.000001)

	ldc, rdc := utils.BinarySearch(func(dc float64) bool {
		parameters.CathodeFallLength = dc
		delta := sourceIntegralAnalyticF(&parameters) - sourceIntegralAnalyticPhelpsF(&parameters) //gammaAnalyticF(parameters)
		return delta > 0
	}, minDc, parameters.SimulationLength(), parameters.CathodeFallLengthPrecision*0.000001)
	return 0.5 * (ldc + rdc)
}

func VoltageRelativeError(sourceIntegral, sourceIntegralMargin, dc, dcMargin float64) float64 {
	return sourceIntegralMargin/(sourceIntegral*(sourceIntegral+1)) + 2.5*(dcMargin/dc)
}

func VoltageCalculation(parameters config.ModelParameters, outputDir, modelName string, logger messages.Logger) (dataRow, altRow utils.ResultInterface, finalModel, altModel *model.Model) {
	minV, maxV := max(100., VoltageFromDc(parameters, parameters.CathodeFallLengthPrecision)), min(1500., VoltageFromDc(parameters, parameters.SimulationLength()-parameters.CathodeFallLengthPrecision))
	stepSize := 10.
	numberOfSteps := int((maxV - minV) / stepSize)
	return GeneralizedCalculation(parameters, numberOfSteps, minV, stepSize,
		func(V float64, mp *config.ModelParameters) {
			mp.CathodeFallPotential = V
			mp.CathodeFallLength = DcFromVoltage(*mp, V)
		},
		func(optResult *OptimizationResult, variance float64, m *model.Model) utils.ResultInterface {
			return &VoltageDataRow{
				OptimizationResult: *optResult,
				VoltageMargin:      m.Parameters.CathodeFallPotential * VoltageRelativeError(optResult.SourceIntegralMonteCarlo, optResult.SourceIntegralMargin, optResult.CathodeFallLength, optResult.CathodeFallLengthMargin),
			}
		}, func(mp *config.ModelParameters, optR *OptimizationResult) bool {
			return mp.GapLength-mp.CathodeFallLength < mp.CathodeFallLengthPrecision || (optR.EffectiveGammaAnalytic > 0 && optR.EffectiveGammaAnalytic > 3*optR.EffectiveGammaMonteCarlo)
		}, false, outputDir, modelName, logger)

}

func VoltageCalculation2(parameters *config.ModelParameters, outputDir, modelName string, logger messages.Logger) (dataRow, altRow utils.ResultInterface, finalModel, altModel *model.Model) {
	minDc, maxDc := DcFromVoltage(*parameters, 80), DcFromVoltage(*parameters, 1500)
	numberOfSteps := int((maxDc - minDc) / parameters.CathodeFallLengthPrecision)
	return GeneralizedCalculation(*parameters, numberOfSteps, minDc, parameters.CathodeFallLengthPrecision,
		func(dc float64, mp *config.ModelParameters) {
			mp.CathodeFallLength = dc
			mp.CathodeFallPotential = VoltageFromDc(*mp, dc)
		},
		func(optResult *OptimizationResult, variance float64, m *model.Model) utils.ResultInterface {
			return &VoltageDataRow{
				OptimizationResult: *optResult,
				VoltageMargin:      m.Parameters.CathodeFallPotential * VoltageRelativeError(optResult.SourceIntegralMonteCarlo, math.Abs(optResult.SourceIntegralDifference)+optResult.SourceIntegralMargin, optResult.CathodeFallLength, optResult.CathodeFallLengthMargin),
			}
		}, func(mp *config.ModelParameters, optR *OptimizationResult) bool {
			return mp.GapLength-mp.CathodeFallLength < mp.CathodeFallLengthPrecision || (optR.EffectiveGammaAnalytic > 0 && optR.EffectiveGammaAnalytic > 3*optR.EffectiveGammaMonteCarlo)
		}, false, outputDir, modelName, logger)

}

type VoltageDataRow struct {
	OptimizationResult
	VoltageMargin float64 `csv:"Voltage margin" units:"Voltage:1"`
}
