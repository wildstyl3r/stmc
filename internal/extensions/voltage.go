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

// var SimplifiedIonDriftVelocityK = map[string][]float64{
// 	"Ar": []float64{47.809144373375744, 25800.},
// 	"He": []float64{240., 79200},
// }

func VoltageFromDc(currentDensity, secondaryEmission, cathodeFallLength, gasDensity float64, species string) float64 {
	lV, rV := utils.BinarySearch(func(V float64) bool {
		return secondaryEmission-gammaAnalyticF(cathodeFallLength, currentDensity, V, gasDensity, utils.IonDriftVelocity[species]) > 0.
	}, 100, 1000, 0.01)
	return 0.5 * (lV + rV)
}

func DcFromVoltage(currentDensity, secondaryEmission, V, gasDensity, precision, gapLength float64, species string) float64 {
	ldc, rdc := utils.BinarySearch(func(dc float64) bool {
		return gammaAnalyticF(dc, currentDensity, V, gasDensity, utils.IonDriftVelocity[species])-secondaryEmission > 0.
	}, precision, gapLength, 0.000001)
	return 0.5 * (ldc + rdc)
}

func voltageCalculationStep(itp *int, dc float64, parameters config.ModelParameters) (loss, gammaIntegral, gammaConfInterval float64, stepModel *model.Model) {
	if itp != nil {
		fmt.Printf("step %d\n", *itp)
		*itp += 1
	} else {
		fmt.Println("preliminary step")
	}
	parameters.CathodeFallLength = dc
	parameters.CathodeFallPotential = VoltageFromDc(
		parameters.CathodeCurrentDensity,
		parameters.SecondaryEmissionCoefficient,
		dc,
		parameters.GasDensity,
		parameters.Species)
	stepModel = model.NewModel(parameters)
	LoadExtensions(stepModel.DataHub)
	stepModel.Run()
	gammaIntegral, gammaConfInterval = gammaIntegralF(stepModel)
	loss = gammaIntegral - parameters.SecondaryEmissionCoefficient
	if parameters.Verbose() {
		fmt.Printf("d_c: %v\nsecondary emission coefficient\n\t integral: %6f\n", dc, gammaIntegral)
	}
	return
}

func VoltageCalculation(configFlags config.Flags, parameters config.ModelParameters, outputDir string) (dataRow *VoltageDataRow, finalModel *model.Model) {
	if *configFlags.Verbose {
		fmt.Print("calculating dc iteratively:\n")
	}
	iteration := 1
	itp := &iteration

	dcMaxLimit := parameters.GapLength
	if parameters.ParallelPlaneHollowCathode {
		dcMaxLimit = parameters.GapLength / 2
	}

	minDc := DcFromVoltage(
		parameters.CathodeCurrentDensity, parameters.SecondaryEmissionCoefficient,
		100,
		parameters.GasDensity, parameters.CathodeFallLengthPrecision, dcMaxLimit, parameters.Species)
	maxDc := min(parameters.GapLength, DcFromVoltage(
		parameters.CathodeCurrentDensity, parameters.SecondaryEmissionCoefficient,
		1000,
		parameters.GasDensity, parameters.CathodeFallLengthPrecision, dcMaxLimit, parameters.Species))

	if *configFlags.Debug {
		return slices.MaxFunc(
			bruteForceStepVoltageCalculation(minDc, maxDc, itp, parameters, outputDir), func(a, b *VoltageDataRow) int {
				return cmp.Compare(math.Abs(a.GammaLoss), math.Abs(b.GammaLoss))
			}), nil
	} else {
		return advancedVoltageCalculation(minDc, maxDc, itp, parameters)
	}
}

type VoltageDataRow struct {
	OptimizationResult
	Voltage       float64 `csv:"Voltage"`
	VoltageMargin float64 `csv:"Voltage"`
}

func bruteForceStepVoltageCalculation(minDc, maxDc float64, itp *int, parameters config.ModelParameters, outputDir string) (debugData []*VoltageDataRow) {
	dcStep := 0.00001
	nSteps := int((maxDc - minDc) / dcStep)
	fmt.Printf("GapLen: %f, nSteps: %d\n", parameters.GapLength, nSteps)
	gamma := make([]float64, nSteps)
	intS := make([]float64, nSteps)
	//debugData := [][]string{{"dc", "E/N", "integrated secondary emission coefficient", "analytic secondary emission coefficient", "gamma difference", "integral gamma CI approximation"}}
	for i := range gamma {
		dc := float64(i)*dcStep + minDc
		var gammaLoss float64
		gammaLoss, gamma[i], _, _ = voltageCalculationStep(itp, dc, parameters)
		intS[i] = 1. / gamma[i]
		V := VoltageFromDc(parameters.CathodeCurrentDensity, parameters.SecondaryEmissionCoefficient, dc, parameters.GasDensity, parameters.Species)
		debugData = append(debugData, &VoltageDataRow{
			OptimizationResult: OptimizationResult{ReducedFieldAtCathode: 2 * V / (dc * parameters.GasDensity * constants.Townsend),
				ReducedFieldAtSheathCenter: V / (dc * parameters.GasDensity * constants.Townsend),
				CathodeFallLength:          dc,
				GammaLoss:                  gammaLoss,
			},
			Voltage: V,
		})
	}
	gammaMean, gammaVariance := utils.MeanAndVariance(gamma, true)
	intSMean, intSVariance := utils.MeanAndVariance(intS, true)
	fmt.Printf("gamma mean: %.9f, gamma variance: %.9f, integral S mean: %.9f, integral S variance: %.9f\n", gammaMean, gammaVariance, intSMean, intSVariance)
	err := utils.WriteAsCSV(debugData, outputDir, "debugGamma")

	if err != nil {
		println("unable to save dc and secondary emission coefficient: ", err.Error())
		os.Exit(1)
	}
	return
}

func advancedVoltageCalculation(minDc, maxDc float64, itp *int, parameters config.ModelParameters) (dataRow *VoltageDataRow, finalModel *model.Model) {
	var gammaLoss, gammaIntegral float64
	gammaI := []float64{}
	var dc float64
	fmt.Println("calculating gamma with stochastic approximation")
	var fLeftLoss, fRightLoss float64
	// var gammaMinDc, gammaMaxDc float64
	{
		fLeftLoss, _, _, _ = voltageCalculationStep(nil, minDc, parameters)
		fRightLoss, _, _, _ = voltageCalculationStep(nil, maxDc, parameters)
	}
	initialDc := 0.5 * (minDc + maxDc) //(parameters.SecondaryEmissionCoefficient - gammaMinDc)/(gammaMaxDc-gammaMinDc)*(maxDc-minDc) + minDc

	approxLossDerivative := (fRightLoss - fLeftLoss) / (maxDc - minDc - 2*parameters.CathodeFallLengthPrecision)

	// var dcConfInterval float64
	minSteps, maxSteps := 10, 700
	dc, _ = utils.StochasticApproximation(minDc, maxDc, initialDc, approxLossDerivative, 0.01, 0.95, true, minSteps, maxSteps, func(dc float64) float64 {
		gammaLoss, gammaIntegral, _, finalModel = voltageCalculationStep(itp, dc, parameters)
		gammaI = append(gammaI, gammaIntegral)
		return gammaLoss
	})
	gammaIntegral = utils.Mean(gammaI[max(len(gammaI)-minSteps, len(gammaI)/4):])
	gammaLoss = parameters.SecondaryEmissionCoefficient - gammaIntegral

	println("saved d_c and γ")
	V := VoltageFromDc(parameters.CathodeCurrentDensity, parameters.SecondaryEmissionCoefficient, dc, parameters.GasDensity, parameters.Species)
	return &VoltageDataRow{
		OptimizationResult: OptimizationResult{ReducedFieldAtCathode: 2 * V / (dc * parameters.GasDensity * constants.Townsend),
			ReducedFieldAtSheathCenter: V / (dc * parameters.GasDensity * constants.Townsend),
			CathodeFallLength:          dc,
			GammaLoss:                  gammaLoss,
		},
		Voltage: V,
	}, finalModel
}
