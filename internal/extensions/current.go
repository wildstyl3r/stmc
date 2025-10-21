package extensions

import (
	"cmp"
	"fmt"
	"math"
	"os"
	"slices"

	"github.com/wildstyl3r/stmc/internal/config"
	"github.com/wildstyl3r/stmc/internal/constants"
	"github.com/wildstyl3r/stmc/internal/model"
	"github.com/wildstyl3r/stmc/internal/utils"
)

func CurrentDensityEquation(dc, secondaryEmissionCoefficient, Vc, N float64, ionDriftVelocity func(float64, float64, float64) float64) float64 {
	return (1 + secondaryEmissionCoefficient) * ionDriftVelocity(Vc, dc, N) * 2 * Vc / (dc * dc) * constants.FreeSpacePermittivityE0
}

func CurrentDensityUncertaintyPropagation(dc, gammaMargin, Vc, N float64, ionDriftVelocity func(float64, float64, float64) float64) float64 {
	return gammaMargin * ionDriftVelocity(Vc, dc, N) * 2 * Vc / (dc * dc) * constants.FreeSpacePermittivityE0
}

func CurrentDensityCalculation(configFlags config.Flags, parameters config.ModelParameters, outputDir string) (dataRow *CurrentDensityDataRow, finalModel *model.Model) {
	if *configFlags.Verbose {
		fmt.Print("calculating dc iteratively:\n")
	}
	iteration := 1
	itp := &iteration
	minDc := parameters.CathodeFallLengthPrecision
	var maxDc float64
	if parameters.ParallelPlaneHollowCathode {
		maxDc = parameters.GapLength / 2
	} else {
		maxDc = parameters.GapLength
	}

	if *configFlags.Debug {
		return slices.MaxFunc(
			bruteForceStepCurrentDensityCalculation(minDc, maxDc, itp, parameters, outputDir), func(a, b *CurrentDensityDataRow) int {
				return cmp.Compare(math.Abs(a.GammaLoss), math.Abs(b.GammaLoss))
			}), nil
	} else {
		return advancedCurrentDensityCalculation(minDc, maxDc, itp, parameters)
	}
}

type CurrentDensityDataRow struct {
	OptimizationResult
	CathodeCurrentDensity float64 `csv:"j Current density"`
	CurrentDensityMargin  float64 `csv:"Current density margin"`
}

func bruteForceStepCurrentDensityCalculation(minDc, maxDc float64, itp *int, parameters config.ModelParameters, outputDir string) (debugData []*CurrentDensityDataRow) {
	dcStep := 0.00001
	nSteps := int((maxDc - minDc) / dcStep)
	fmt.Printf("GapLen: %f, nSteps: %d\n", parameters.GapLength, nSteps)
	gamma := make([]float64, nSteps)
	intS := make([]float64, nSteps)
	for i := range gamma {
		dc := float64(i)*dcStep + minDc
		var gammaLoss, gammaMargin float64
		gammaLoss, gamma[i], gammaMargin, _ = currentDensityCalculationStep(itp, dc, parameters)
		currentDensity := CurrentDensityEquation(dc, parameters.SecondaryEmissionCoefficient, parameters.CathodeFallPotential, parameters.GasDensity, utils.IonDriftVelocity[parameters.Species])
		currentDensityMargin := CurrentDensityUncertaintyPropagation(dc, gammaMargin, parameters.CathodeFallPotential, parameters.GasDensity, utils.IonDriftVelocity[parameters.Species])
		intS[i] = 1. / gamma[i]
		debugData = append(debugData, &CurrentDensityDataRow{
			OptimizationResult: OptimizationResult{ReducedFieldAtCathode: 2 * parameters.CathodeFallPotential / (dc * parameters.GasDensity * constants.Townsend),
				ReducedFieldAtSheathCenter: parameters.CathodeFallPotential / (dc * parameters.GasDensity * constants.Townsend),
				CathodeFallLength:          dc,
				GammaLoss:                  gammaLoss,
			},
			CathodeCurrentDensity: currentDensity,
			CurrentDensityMargin:  currentDensityMargin,
		})
	}
	gammaMean, gammaVariance := utils.MeanAndVariance(gamma, true)
	intSMean, intSVariance := utils.MeanAndVariance(intS, true)
	fmt.Printf("gamma mean: %.9f, gamma variance: %.9f, integral S mean: %.9f, integral S variance: %.9f\n", gammaMean, gammaVariance, intSMean, intSVariance)

	err := utils.WriteAsCSV(debugData, outputDir, "debugCurrent")

	if err != nil {
		println("unable to save dc and current: ", err.Error())
		os.Exit(1)
	}
	return
}

func advancedCurrentDensityCalculation(minDc, maxDc float64, itp *int, parameters config.ModelParameters) (dataRow *CurrentDensityDataRow, finalModel *model.Model) {
	var gammaLoss, gammaIntegral, gammaMargin float64
	gammaI := []float64{}
	var dc float64
	fmt.Println("calculating current density with stochastic approximation")
	var fLeftLoss, fRightLoss float64
	{
		fLeftLoss, _, _, _ = currentDensityCalculationStep(nil, minDc+parameters.CathodeFallLengthPrecision, parameters)
		fRightLoss, _, _, _ = currentDensityCalculationStep(nil, maxDc-parameters.CathodeFallLengthPrecision, parameters)
	}
	approxLossDerivative := (fRightLoss - fLeftLoss) / (maxDc - minDc - 2*parameters.CathodeFallLengthPrecision)

	minSteps := 10
	maxSteps := 700
	dc, _ = utils.StochasticApproximation(minDc, maxDc, 0.5*(maxDc+minDc), approxLossDerivative, 0.002, 0.95, true, minSteps, maxSteps, func(dc float64) float64 {
		gammaLoss, gammaIntegral, gammaMargin, finalModel = currentDensityCalculationStep(itp, dc, parameters)
		gammaI = append(gammaI, gammaIntegral)
		return gammaLoss
	})
	gammaIntegral = utils.Mean(gammaI[max(len(gammaI)-minSteps, len(gammaI)/4):])
	gammaMargin = max(gammaMargin, utils.StudentedMarginFromData(0.95, gammaI[max(len(gammaI)-minSteps, len(gammaI)/4):]))
	gammaLoss = parameters.SecondaryEmissionCoefficient - gammaIntegral

	currentDensity := CurrentDensityEquation(dc, parameters.SecondaryEmissionCoefficient, parameters.CathodeFallPotential, parameters.GasDensity, utils.IonDriftVelocity[parameters.Species])
	currentDensityMargin := CurrentDensityUncertaintyPropagation(dc, gammaMargin, parameters.CathodeFallPotential, parameters.GasDensity, utils.IonDriftVelocity[parameters.Species])
	finalModel.Parameters.CathodeCurrentDensity = currentDensity

	println("saved d_c and j")
	return &CurrentDensityDataRow{
		OptimizationResult: OptimizationResult{ReducedFieldAtCathode: 2 * parameters.CathodeFallPotential / (dc * parameters.GasDensity * constants.Townsend),
			ReducedFieldAtSheathCenter: parameters.CathodeFallPotential / (dc * parameters.GasDensity * constants.Townsend),
			CathodeFallLength:          dc,
			GammaLoss:                  gammaLoss,
		},
		CathodeCurrentDensity: currentDensity,
		CurrentDensityMargin:  currentDensityMargin,
	}, finalModel
}

func currentDensityCalculationStep(itp *int, dc float64, parameters config.ModelParameters) (loss, gammaIntegral, gammaConfInterval float64, stepModel *model.Model) {
	if itp != nil {
		fmt.Printf("step %d\n", *itp)
		*itp += 1
	} else {
		fmt.Println("preliminary step")
	}
	parameters.CathodeFallLength = dc
	stepModel = (model.NewModel(parameters))
	LoadExtensions(stepModel.DataHub)
	gammaIntegral, gammaConfInterval = gammaIntegralF(stepModel, 0)
	loss = parameters.SecondaryEmissionCoefficient - gammaIntegral
	if parameters.Verbose() {
		fmt.Printf("d_c: %v\nsecondary emission coefficient\n\t integral: %6f\n\t target:%6f\n", dc, gammaIntegral, parameters.SecondaryEmissionCoefficient)
	}
	return
}
