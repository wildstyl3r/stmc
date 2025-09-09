package extensions

import (
	"encoding/csv"
	"fmt"
	"math"
	"os"
	"strconv"

	"github.com/wildstyl3r/lxgata"
	"github.com/wildstyl3r/stmc/internal/config"
	"github.com/wildstyl3r/stmc/internal/constants"
	"github.com/wildstyl3r/stmc/internal/datahub"
	"github.com/wildstyl3r/stmc/internal/model"
	"github.com/wildstyl3r/stmc/internal/utils"
)

func gammaIntegralF(model *model.Model) float64 {
	ionizations := utils.GriddedInterval{
		Step:   model.XStep,
		Values: datahub.Get(datahub.KeyType(lxgata.IONIZATION), model).([]float64),
	}
	sourceTermIntegral, err := utils.TrapezoidIntegration(0, float64(utils.Argmax(datahub.Get("GlowDischargeDensity", model).([]float64)))*ionizations.Step, ionizations)
	if err != nil {
		fmt.Printf("error calculating integral gamma: %#v", err)
	}
	sourceTermIntegral *= model.Parameters.Pressure //* Torr / model.Parameters.Pressure / de.cathodeFlux
	return 1 / sourceTermIntegral

}

// not applicable for UniformField
func gammaAnalyticF(dc, j, Vc, N float64, ionDriftVelocity func(float64, float64, float64) float64) float64 {
	return j*dc*dc/(2.*Vc*ionDriftVelocity(Vc, dc, N)*constants.FreeSpacePermittivityE0) - 1. //-> 0
}

func gammaWithLoss(model *model.Model, loss utils.LossType) (lossValue float64, gi float64, ga float64) {
	ga = gammaAnalyticF(
		model.Parameters.CathodeFallLength,
		model.Parameters.CathodeCurrentDensity,
		model.Parameters.CathodeFallPotential,
		model.Parameters.GasDensity,
		utils.IonDriftVelocity[model.Parameters.Species])
	gi = gammaIntegralF(model)
	switch loss {
	case utils.MSE:
		lossValue = (ga - gi) * (ga - gi)
	case utils.Difference:
		lossValue = ga - gi
	}
	return
}

func EstimateCathodeFallLengthLimits(parameters *config.ModelParameters) (from float64, to float64) {
	from, _ = utils.BinarySearch(func(dc float64) bool {
		return 0. < gammaAnalyticF(dc,
			parameters.CathodeCurrentDensity,
			parameters.CathodeFallPotential,
			parameters.GasDensity,
			utils.IonDriftVelocity[parameters.Species])
	}, 1e-8, parameters.GapLength, parameters.CathodeFallLengthPrecision*0.01)

	_, to = utils.BinarySearch(func(dc float64) bool {
		return 1 < gammaAnalyticF(dc,
			parameters.CathodeCurrentDensity,
			parameters.CathodeFallPotential,
			parameters.GasDensity,
			utils.IonDriftVelocity[parameters.Species]) // might be false everywhere in the gap, but as close as possible to true domain
	}, 0, parameters.GapLength, parameters.CathodeFallLengthPrecision*0.01)
	return from, to
}

func getApproximateCathodeFallLengthForGamma(gamma, minDc, maxDc float64, parameters config.ModelParameters) float64 {
	initialDcL, initialDcR := utils.BinarySearch(func(dc float64) bool {
		return gamma < gammaAnalyticF(
			dc,
			parameters.CathodeCurrentDensity,
			parameters.CathodeFallPotential,
			parameters.GasDensity,
			utils.IonDriftVelocity[parameters.Species])
	}, minDc, maxDc, parameters.CathodeFallLengthPrecision*0.1)
	return 0.5 * (initialDcL + initialDcR)
}

func gammaCalculationStep(itp *int, dc float64, parameters config.ModelParameters, lossType utils.LossType) (loss, gammaIntegral, gammaAnalytic float64) {
	if itp != nil {
		fmt.Printf("step %d\n", *itp)
		*itp += 1
	} else {
		fmt.Println("preliminary step")
	}
	parameters.GapLength = dc
	model := model.NewModel(parameters)
	model.Run()
	datahub.Reset(&model)
	loss, gammaIntegral, gammaAnalytic = gammaWithLoss(&model, lossType)
	if parameters.Verbose() {
		fmt.Printf("d_c: %v\nsecondary emission coefficient\n\t integral: %6f\n\t analytic:%6f\n", dc, gammaIntegral, gammaAnalytic)
	}
	return
}

func GammaCalculation(configFlags config.Flags, parameters config.ModelParameters, modelName string, outputDir string) []string {
	if *configFlags.Verbose {
		fmt.Print("calculating dc iteratively:\n")
	}
	iteration := 1
	itp := &iteration
	minDc, maxDc := EstimateCathodeFallLengthLimits(&parameters)
	if *configFlags.Debug {
		return bruteForceStepGammaCalculation(minDc, maxDc, configFlags, itp, parameters, modelName, outputDir)
	} else {
		return advancedGammaCalculation(minDc, maxDc, configFlags, itp, parameters)
	}
}

func bruteForceStepGammaCalculation(minDc, maxDc float64, configFlags config.Flags, itp *int, parameters config.ModelParameters, modelName, outputDir string) []string {
	dcStep := 0.00001
	nSteps := int((maxDc - minDc) / dcStep)
	fmt.Printf("GapLen: %f, nSteps: %d\n", parameters.GapLength, nSteps)
	gamma := make([]float64, nSteps)
	intS := make([]float64, nSteps)
	debugData := [][]string{{"dc", "E/N", "integrated secondary emission coefficient", "analytic secondary emission coefficient", "gamma difference"}}
	for i := range gamma {
		dc := float64(i)*dcStep + minDc
		var lossType utils.LossType
		switch *configFlags.RootFindingAlgorithm {
		case "t":
			lossType = utils.MSE
		case "b", "s":
			lossType = utils.Difference
		}
		var gammaAnalytic, gammaLoss float64
		gammaLoss, gamma[i], gammaAnalytic = gammaCalculationStep(itp, dc, parameters, lossType)
		intS[i] = 1. / gamma[i]
		debugData = append(debugData, []string{
			strconv.FormatFloat(dc, 'f', 10, 64),
			strconv.FormatFloat(parameters.CathodeFallPotential/(dc*parameters.GasDensity*constants.Townsend), 'f', 10, 64),
			strconv.FormatFloat(gamma[i], 'f', 10, 64),
			strconv.FormatFloat(gammaAnalytic, 'f', 10, 64),
			strconv.FormatFloat(gammaLoss, 'f', 10, 64),
		})
	}
	gammaMean, gammaVariance := utils.MeanAndVariance(gamma, true)
	intSMean, intSVariance := utils.MeanAndVariance(intS, true)
	fmt.Printf("[%s] gamma mean: %.9f, gamma variance: %.9f, integral S mean: %.9f, integral S variance: %.9f\n", modelName, gammaMean, gammaVariance, intSMean, intSVariance)
	debugFile, err := utils.OpenFile(true, outputDir, "debug", utils.GetFilename(*configFlags.ConfigFileNamePointer)+modelName, utils.TypeTXT)

	if err != nil {
		println("unable to save dc and secondary emission coefficient: ", err.Error())
		os.Exit(1)
	} else {
		debugW := csv.NewWriter(debugFile)
		debugW.WriteAll(debugData)
		debugW.Flush()
	}
	return nil
}

func advancedGammaCalculation(minDc, maxDc float64, configFlags config.Flags, itp *int, parameters config.ModelParameters) []string {
	var gammaLoss, gammaIntegral, gammaAnalytic, gammaCI float64
	var dc float64
	switch *configFlags.RootFindingAlgorithm {
	case "s":
		fmt.Println("calculating gamma with stochastic approximation")
		var fLeftLoss, fRightLoss, initialDc float64
		{
			var gILeft, gIRight float64
			fLeftLoss, gILeft, _ = gammaCalculationStep(nil, minDc+parameters.CathodeFallLengthPrecision, parameters, utils.Difference)
			fRightLoss, gIRight, _ = gammaCalculationStep(nil, maxDc-parameters.CathodeFallLengthPrecision, parameters, utils.Difference)

			meanGI := 0.5 * (gILeft + gIRight)
			initialDc = getApproximateCathodeFallLengthForGamma(meanGI, minDc, maxDc, parameters)
		}
		approxLossDerivative := (fRightLoss - fLeftLoss) / (maxDc - minDc - 2*parameters.CathodeFallLengthPrecision)

		gammaI := []float64{}
		dc = utils.StochasticApproximation(minDc, maxDc, initialDc, approxLossDerivative, parameters.CathodeFallLengthPrecision, constants.Quantile95, 10, func(dc float64) float64 {
			gammaLoss, gammaIntegral, gammaAnalytic = gammaCalculationStep(itp, dc, parameters, utils.Difference)
			gammaI = append(gammaI, gammaIntegral)
			return gammaLoss
		})
		var gammaVariance float64
		gammaIntegral, gammaVariance = utils.MeanAndVariance(gammaI[len(gammaI)-10:], true)
		gammaCI = 2 * math.Sqrt(gammaVariance) * constants.Quantile95 / math.Sqrt(float64(10))
		gammaLoss = gammaAnalytic - gammaIntegral

	case "b":
		fmt.Println("calculating gamma with naive bisection")
		dcLeft, dcRight := utils.BinarySearch(func(dc float64) bool {
			gammaLoss, gammaIntegral, gammaAnalytic = gammaCalculationStep(itp, dc, parameters, utils.Difference)
			return gammaLoss > 0
		}, minDc, min(maxDc, parameters.GapLength), parameters.CathodeFallLengthPrecision)
		dc = 0.5 * (dcLeft + dcRight)

	case "t":
		fmt.Println("calculating gamma with naive ternary search")
		fmt.Printf("min dc: %f, max dc: %f, gap len: %f\n", minDc, maxDc, parameters.GapLength)
		dc = utils.TernarySearchMax(func(dc float64) float64 {
			gammaLoss, gammaIntegral, gammaAnalytic = gammaCalculationStep(itp, dc, parameters, utils.MSE)
			return -gammaLoss
		}, minDc, min(maxDc, parameters.GapLength), parameters.CathodeFallLengthPrecision)
	}

	println("saved d_c and γ")
	return []string{
		strconv.FormatFloat(parameters.CathodeFallPotential/(dc*parameters.GasDensity*constants.Townsend), 'f', 10, 64),
		strconv.FormatFloat(gammaIntegral, 'f', 10, 64),
		strconv.FormatFloat(gammaAnalytic, 'f', 10, 64),
		strconv.FormatFloat(gammaLoss, 'f', 10, 64),
		strconv.FormatFloat(dc, 'f', 10, 64),
		strconv.FormatFloat(gammaCI, 'f', 10, 64),
	}
}
