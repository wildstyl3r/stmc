package extensions

import (
	"fmt"
	"math"
	"os"

	"github.com/wildstyl3r/stmc/internal/config"
	"github.com/wildstyl3r/stmc/internal/constants"
	"github.com/wildstyl3r/stmc/internal/model"
	"github.com/wildstyl3r/stmc/internal/utils"
)

func sourceIntegralMonteCarloF(m *model.Model, precision float64) (sumMean, sumVariance float64) {
	m.Run(func(m *model.Model) int {
		if m.TotalElectronsPassed == 0 {
			return m.Parameters.NElectrons
		} else {
			sumMean, sumVariance = m.GetMetrics(SourceIntegralPerCathodeElectronFluxKey).(float64), m.GetMetrics(SourceIntegralPerCathodeElectronFluxVarianceKey).(float64)
			sumMargin := utils.StudentedMargin(0.95, sumVariance, m.TotalElectronsPassed)
			if precision < sumMargin*2 && precision != 0 {
				return m.Parameters.AddByNElectrons
			} else {
				fmt.Printf("one step ion source to max(n) integral margin: %f\n", sumMargin)
				return 0
			}
		}
	})
	return sumMean, sumVariance
}

// not applicable for UniformField
func sourceIntegralAnalyticF(dc, j, Vc, N float64, ionDriftVelocity func(float64, float64, float64) float64) float64 {
	return 1. / (j*dc*dc/(2.*Vc*ionDriftVelocity(Vc, dc, N)*constants.FreeSpacePermittivityE0) - 1.) //-> 0
}

func sourceWithLoss(model *model.Model, sourcePrecision float64) (lossValue, Si, Sa, SiVariance float64) {
	Sa = sourceIntegralAnalyticF(
		model.Parameters.CathodeFallLength,
		model.Parameters.CathodeCurrentDensity,
		model.Parameters.CathodeFallPotential,
		model.Parameters.GasDensity,
		utils.IonDriftVelocity[model.Parameters.Species])
	Si, SiVariance = sourceIntegralMonteCarloF(model, sourcePrecision)
	lossValue = Sa - Si
	return
}

func getApproximateCathodeFallLengthForSourceIntegral(sourceIntegral, minDc, maxDc float64, parameters config.ModelParameters) float64 {
	initialDcL, initialDcR := utils.BinarySearch(func(dc float64) bool {
		si := sourceIntegralAnalyticF(
			dc,
			parameters.CathodeCurrentDensity,
			parameters.CathodeFallPotential,
			parameters.GasDensity,
			utils.IonDriftVelocity[parameters.Species])
		return sourceIntegral > si
	}, minDc, maxDc, parameters.CathodeFallLengthPrecision*0.001)
	return 0.5 * (initialDcL + initialDcR)
}

func sourceIntegralCalculationStep(itp, itn *int, dc, sourcePrecision float64, parameters config.ModelParameters) (loss, sourceIntegralMonteCarlo, sourceIntegralAnalytic, sourceIntegralVariance float64, stepModel *model.Model) {
	if itp != nil {
		if itn != nil {
			fmt.Printf("step %d of %d\n", *itp, *itn)
		} else {
			fmt.Printf("step %d\n", *itp)
		}
		*itp += 1
	} else {
		fmt.Println("preliminary step")
	}
	parameters.CathodeFallLength = dc
	stepModel = model.NewModel(parameters)
	LoadExtensions(stepModel.DataHub)
	loss, sourceIntegralMonteCarlo, sourceIntegralAnalytic, sourceIntegralVariance = sourceWithLoss(stepModel, sourcePrecision)
	if parameters.Verbose() {
		fmt.Printf("d_c: %v\nion source integral \n\t Monte Carlo: %6f\n\t analytic:%6f\n", dc, sourceIntegralMonteCarlo, sourceIntegralAnalytic)
	}
	return
}

func SourceIntegralCalculation(configFlags config.Flags, parameters config.ModelParameters, outputDir, modelName string) (dataRow *SourceIntegralDataRow, finalModel *model.Model) {
	if *configFlags.Verbose {
		fmt.Print("calculating dc iteratively:\n")
	}
	iteration := 1
	itp := &iteration
	minDc, maxDc := EstimateCathodeFallLengthLimits(&parameters)
	if *configFlags.Debug {
		return bruteForceStepSourceIntegralCalculation(minDc, maxDc, itp, parameters, outputDir, modelName)
	} else {
		return advancedSourceIntegralCalculation(minDc, maxDc, itp, parameters)
	}
}

type SourceIntegralDataRow struct {
	OptimizationResult
	SourceIntegralAnalytic   float64 `csv:"Source integral analytic"`
	SourceIntegralMonteCarlo float64 `csv:"Source integral Monte Carlo"`
	SourceIntegralLoss       float64 `csv:"Source integral loss"`
	SourceIntegralMargin     float64 `csv:"Source integral margin"`
	GammaIntegral            float64 `csv:"Gamma integral"`
	GammaIntegralMargin      float64 `csv:"Gamma integral margin"`
	GammaAnalytic            float64 `csv:"Gamma analytic"`
}

func bruteForceStepSourceIntegralCalculation(minDc, maxDc float64, itp *int, parameters config.ModelParameters, outputDir, modelName string) (dataRow *SourceIntegralDataRow, finalModel *model.Model) {
	dcStep := 0.00001
	nSteps := int((maxDc - minDc) / dcStep)
	fmt.Printf("GapLen: %f, nSteps: %d\n", parameters.GapLength, nSteps)
	sourceIntegral := make([]float64, nSteps)
	intS := make([]float64, nSteps)
	//debugData := [][]string{{"dc", "E/N", "integrated secondary emission coefficient", "analytic secondary emission coefficient", "gamma difference", "integral gamma CI approximation"}}
	bestLoss := math.Inf(+1)
	var debugData []*SourceIntegralDataRow
	for i := range sourceIntegral {
		dc := float64(i)*dcStep + minDc
		var sourceIntegralAnalytic, sourceIntegralLoss, sourceIntegralMonteCarloMargin float64
		var m *model.Model
		sourceIntegralLoss, sourceIntegral[i], sourceIntegralAnalytic, sourceIntegralMonteCarloMargin, m = sourceIntegralCalculationStep(itp, &nSteps, dc, 0, parameters)

		intS[i] = 1. / sourceIntegral[i]
		if sourceIntegralAnalytic > 0 {
			debugData = append(debugData, &SourceIntegralDataRow{
				OptimizationResult: OptimizationResult{ReducedFieldAtCathode: m.ReducedFieldAtCathode(),
					ReducedFieldAtSheathCenter: m.ReducedFieldMidSheath(),
					CathodeFallLength:          dc,
					GammaLoss:                  -sourceIntegralLoss / (sourceIntegral[i] * sourceIntegralAnalytic),
				},
				SourceIntegralAnalytic:   sourceIntegralAnalytic,
				SourceIntegralMonteCarlo: sourceIntegral[i],
				SourceIntegralLoss:       sourceIntegralLoss,
				SourceIntegralMargin:     sourceIntegralMonteCarloMargin,
			})
			if math.Abs(sourceIntegralLoss) < math.Abs(bestLoss) {
				bestLoss = sourceIntegralLoss
				finalModel = m
				dataRow = debugData[len(debugData)-1]
			}
		}

	}
	// gammaMean, gammaVariance := utils.MeanAndVariance(sourceIntegral, true)
	// intSMean, intSVariance := utils.MeanAndVariance(intS, true)
	// fmt.Printf("gamma mean: %.9f, gamma variance: %.9f, integral S mean: %.9f, integral S variance: %.9f\n", gammaMean, gammaVariance, intSMean, intSVariance)
	err := utils.WriteAsCSV(debugData, outputDir+"/"+modelName, "source_integral_bruteforce")

	if err != nil {
		println("unable to save dc and secondary emission coefficient: ", err.Error())
		os.Exit(1)
	}
	return
}

func advancedSourceIntegralCalculation(minDc, maxDc float64, itp *int, parameters config.ModelParameters) (dataRow *SourceIntegralDataRow, finalModel *model.Model) {
	var sourceIntegralLoss, sourceIntegralMonteCarlo, sourceIntegralAnalytic float64
	sourceIntegralMC := []float64{}
	sourceIntegralMCVariance := []float64{}
	var dc float64
	fmt.Println("calculating source integral up to x_{max(n)} with stochastic approximation")
	var fLeftLoss, fRightLoss, initialDc float64
	{
		var sIMCLeft, sIMCRight float64
		fLeftLoss, sIMCLeft, _, _, _ = sourceIntegralCalculationStep(nil, nil, minDc, 0.7, parameters)
		fRightLoss, sIMCRight, _, _, _ = sourceIntegralCalculationStep(nil, nil, maxDc, 0.7, parameters)

		k := (sIMCRight - sIMCLeft) / (maxDc - minDc)
		b := sIMCRight - k*maxDc

		dcLeft, dcRight := utils.BinarySearch(func(dc float64) bool {
			delta := sourceIntegralAnalyticF(dc, parameters.CathodeCurrentDensity, parameters.CathodeFallPotential, parameters.GasDensity, utils.IonDriftVelocity[parameters.Species]) -
				(k*dc + b)
			return delta < 0
		}, minDc, maxDc, parameters.CathodeFallLengthPrecision)

		initialDc = 0.5 * (dcLeft + dcRight)
	}
	approxLossDerivative := (fRightLoss - fLeftLoss) / (maxDc - minDc)

	// var dcConfInterval float64
	minSteps, maxSteps := 40, 200
	var suggestion float64
	sourceIntegralPrecision := 0.1
	dc, _ = utils.StochasticApproximationWithSameSizeSamples(minDc, maxDc, initialDc, approxLossDerivative, sourceIntegralPrecision, 0.95, true, minSteps, maxSteps, parameters.NElectrons, func(dc float64) (lossMean, lossVariance float64) {
		sourceIntegralLoss, sourceIntegralMonteCarlo, sourceIntegralAnalytic, lossVariance, finalModel = sourceIntegralCalculationStep(itp, nil, dc, 0, parameters)
		sourceIntegralMC = append(sourceIntegralMC, sourceIntegralMonteCarlo)
		sourceIntegralMCVariance = append(sourceIntegralMCVariance, lossVariance)
		return sourceIntegralLoss, lossVariance
	}, func() *float64 {
		if *itp%5 == 0 && *itp > 0 { //&& math.Abs(utils.LinearlyWeightedSameSizePoolMean(sourceIntegralMC[len(sourceIntegralMC)-4:])-sourceIntegralAnalytic) > sourceIntegralPrecision*0.25 {
			backsight := 4
			if *itp%10 == 0 {
				backsight = 9
			}
			suggestion = getApproximateCathodeFallLengthForSourceIntegral(utils.Mean(sourceIntegralMC[len(sourceIntegralMC)-backsight:]), minDc, maxDc, parameters) //sourceIntegralMonteCarlo
			return &suggestion
		} else {
			return nil
		}
	})
	sourceIntegralAnalytic = sourceIntegralAnalyticF(dc, parameters.CathodeCurrentDensity, parameters.CathodeFallPotential, parameters.GasDensity, utils.IonDriftVelocity[parameters.Species])
	subsampleLen := max(len(sourceIntegralMC)-minSteps, len(sourceIntegralMC)/3)
	sourceIntegralMonteCarlo = utils.Mean(sourceIntegralMC[subsampleLen:])
	finalMCVariance := utils.SameSizePoolVariance(sourceIntegralMCVariance[subsampleLen:])
	sourceIntegralMargin := utils.StudentedMargin(0.95, finalMCVariance, subsampleLen*parameters.NElectrons)
	sourceIntegralLoss = sourceIntegralAnalytic - sourceIntegralMonteCarlo

	println("saved d_c and γ")
	return &SourceIntegralDataRow{
		OptimizationResult: OptimizationResult{
			ReducedFieldAtCathode:      finalModel.ReducedFieldAtCathode(),
			ReducedFieldAtSheathCenter: finalModel.ReducedFieldMidSheath(),
			CathodeFallLength:          dc,
			GammaLoss:                  sourceIntegralLoss,
		},
		SourceIntegralAnalytic:   sourceIntegralAnalytic,
		SourceIntegralMonteCarlo: sourceIntegralMonteCarlo,
		SourceIntegralLoss:       sourceIntegralLoss,
		SourceIntegralMargin:     sourceIntegralMargin,
		GammaIntegral:            1. / sourceIntegralMonteCarlo,
		GammaIntegralMargin:      (sourceIntegralMargin + math.Abs(sourceIntegralLoss)) / (sourceIntegralMonteCarlo * sourceIntegralMonteCarlo),
		GammaAnalytic:            1. / sourceIntegralAnalytic,
	}, finalModel
}
