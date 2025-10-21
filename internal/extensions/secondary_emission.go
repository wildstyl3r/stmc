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

func gammaIntegralF(m *model.Model, precision float64) (gamma, gammaMargin float64) {
	m.Run(func(m *model.Model) int {
		if m.TotalElectronsPassed == 0 {
			return m.Parameters.NElectrons
		} else {
			sumMean, sumVariance := m.GetMetrics(SourceIntegralPerCathodeElectronFluxKey).(float64), m.GetMetrics(SourceIntegralPerCathodeElectronFluxVarianceKey).(float64)
			sumMargin := utils.StudentedMargin(0.95, sumVariance, m.Parameters.NElectrons)
			gamma = 1 / sumMean
			gammaMargin = sumMargin / (sumMean * sumMean)
			fmt.Printf("one step gamma CI: %f\n", gammaMargin)
			if precision < 2*gammaMargin && precision != 0 {
				return m.Parameters.AddByNElectrons
			} else {
				return 0
			}
		}
	})

	return gamma, gammaMargin
}

// not applicable for UniformField
func gammaAnalyticF(dc, j, Vc, N float64, ionDriftVelocity func(float64, float64, float64) float64) float64 {
	return j*dc*dc/(2.*Vc*ionDriftVelocity(Vc, dc, N)*constants.FreeSpacePermittivityE0) - 1. //-> 0
}

func gammaWithLoss(model *model.Model) (lossValue, gi, ga, giMargin float64) {
	ga = gammaAnalyticF(
		model.Parameters.CathodeFallLength,
		model.Parameters.CathodeCurrentDensity,
		model.Parameters.CathodeFallPotential,
		model.Parameters.GasDensity,
		utils.IonDriftVelocity[model.Parameters.Species])
	gi, giMargin = gammaIntegralF(model, 0)
	lossValue = ga - gi
	return
}

func EstimateCathodeFallLengthLimits(parameters *config.ModelParameters) (from float64, to float64) {
	var upperBound float64
	if parameters.ParallelPlaneHollowCathode {
		upperBound = parameters.GapLength / 2
	} else {
		upperBound = parameters.GapLength
	}
	upperBound -= parameters.CathodeFallLengthPrecision * 0.01

	lowerBound := 1e-9

	from, _ = utils.BinarySearch(func(dc float64) bool {
		g := gammaAnalyticF(dc,
			parameters.CathodeCurrentDensity,
			parameters.CathodeFallPotential,
			parameters.GasDensity,
			utils.IonDriftVelocity[parameters.Species])
		return 0. < g
	}, lowerBound, upperBound, parameters.CathodeFallLengthPrecision*0.01)

	_, to = utils.BinarySearch(func(dc float64) bool {
		g := gammaAnalyticF(dc,
			parameters.CathodeCurrentDensity,
			parameters.CathodeFallPotential,
			parameters.GasDensity,
			utils.IonDriftVelocity[parameters.Species])
		return 1 < g // might be false everywhere in the gap, but as close as possible to true domain
	}, lowerBound, upperBound, parameters.CathodeFallLengthPrecision*0.01)
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

func gammaCalculationStep(itp *int, dc float64, parameters config.ModelParameters) (loss, gammaIntegral, gammaAnalytic, gammaIntegralMargin float64, stepModel *model.Model) {
	if itp != nil {
		fmt.Printf("step %d\n", *itp)
		*itp += 1
	} else {
		fmt.Println("preliminary step")
	}
	parameters.CathodeFallLength = dc
	stepModel = model.NewModel(parameters)
	LoadExtensions(stepModel.DataHub)
	loss, gammaIntegral, gammaAnalytic, gammaIntegralMargin = gammaWithLoss(stepModel)
	if parameters.Verbose() {
		fmt.Printf("d_c: %v\nsecondary emission coefficient\n\t integral: %6f\n\t analytic:%6f\n", dc, gammaIntegral, gammaAnalytic)
	}
	return
}

func GammaCalculation(configFlags config.Flags, parameters config.ModelParameters, outputDir string) (dataRow *GammaDataRow, finalModel *model.Model) {
	if *configFlags.Verbose {
		fmt.Print("calculating dc iteratively:\n")
	}
	iteration := 1
	itp := &iteration
	minDc, maxDc := EstimateCathodeFallLengthLimits(&parameters)
	if *configFlags.Debug {
		return slices.MaxFunc(
			bruteForceStepGammaCalculation(minDc, maxDc, itp, parameters, outputDir), func(a, b *GammaDataRow) int {
				return cmp.Compare(math.Abs(a.GammaLoss), math.Abs(b.GammaLoss))
			}), nil
	} else {
		return advancedGammaCalculation(minDc, maxDc, itp, parameters)
	}
}

type OptimizationResult struct {
	ModelName                  string  `csv:"Model name"`
	ReducedFieldAtCathode      float64 `csv:"E/n at cathode"`
	ReducedFieldAtSheathCenter float64 `csv:"E/n at mid-sheath"`
	CathodeFallLength          float64 `csv:"Sheath length"`
	GammaLoss                  float64 `csv:"Gamma loss"`
}

func (row *OptimizationResult) Index() float64 {
	return row.ReducedFieldAtCathode
}

type GammaDataRow struct {
	OptimizationResult
	GammaAnalytic float64 `csv:"Gamma analytic"`
	GammaIntegral float64 `csv:"Gamma integral"`
	GammaMargin   float64 `csv:"Gamma confidence interval"`
}

func bruteForceStepGammaCalculation(minDc, maxDc float64, itp *int, parameters config.ModelParameters, outputDir string) (debugData []*GammaDataRow) {
	dcStep := 0.00001
	nSteps := int((maxDc - minDc) / dcStep)
	fmt.Printf("GapLen: %f, nSteps: %d\n", parameters.GapLength, nSteps)
	gamma := make([]float64, nSteps)
	intS := make([]float64, nSteps)
	//debugData := [][]string{{"dc", "E/N", "integrated secondary emission coefficient", "analytic secondary emission coefficient", "gamma difference", "integral gamma CI approximation"}}
	for i := range gamma {
		dc := float64(i)*dcStep + minDc
		var gammaAnalytic, gammaLoss, integralGammaMargin float64
		gammaLoss, gamma[i], gammaAnalytic, integralGammaMargin, _ = gammaCalculationStep(itp, dc, parameters)
		intS[i] = 1. / gamma[i]
		debugData = append(debugData, &GammaDataRow{
			OptimizationResult: OptimizationResult{ReducedFieldAtCathode: 2 * parameters.CathodeFallPotential / (dc * parameters.GasDensity * constants.Townsend),
				ReducedFieldAtSheathCenter: parameters.CathodeFallPotential / (dc * parameters.GasDensity * constants.Townsend),
				CathodeFallLength:          dc,
				GammaLoss:                  gammaLoss,
			},
			GammaAnalytic: gammaAnalytic,
			GammaIntegral: gamma[i],
			GammaMargin:   integralGammaMargin,
		})
	}
	gammaMean, gammaVariance := utils.MeanAndVariance(gamma, true)
	intSMean, intSVariance := utils.MeanAndVariance(intS, true)
	fmt.Printf("gamma mean: %.9f, gamma variance: %.9f, integral S mean: %.9f, integral S variance: %.9f\n", gammaMean, gammaVariance, intSMean, intSVariance)
	err := utils.WriteAsCSV(debugData, outputDir, "gamma_bruteforce")

	if err != nil {
		println("unable to save dc and secondary emission coefficient: ", err.Error())
		os.Exit(1)
	}
	return
}

func advancedGammaCalculation(minDc, maxDc float64, itp *int, parameters config.ModelParameters) (dataRow *GammaDataRow, finalModel *model.Model) {
	var gammaLoss, gammaIntegral, gammaAnalytic, gammaMargin, lastStepGammaConfidenceInterval float64
	gammaI := []float64{}
	var dc float64
	fmt.Println("calculating gamma with stochastic approximation")
	var fLeftLoss, fRightLoss, initialDc float64
	{
		var gILeft, gIRight float64
		fLeftLoss, gILeft, _, _, _ = gammaCalculationStep(nil, minDc+parameters.CathodeFallLengthPrecision, parameters)
		fRightLoss, gIRight, _, _, _ = gammaCalculationStep(nil, maxDc-parameters.CathodeFallLengthPrecision, parameters)

		meanGI := 0.5 * (gILeft + gIRight)
		initialDc = getApproximateCathodeFallLengthForGamma(meanGI, minDc, maxDc, parameters)
	}
	approxLossDerivative := (fRightLoss - fLeftLoss) / (maxDc - minDc - 2*parameters.CathodeFallLengthPrecision)

	// var dcConfInterval float64
	minSteps, maxSteps := 10, 700
	dc, _ = utils.StochasticApproximation(minDc, maxDc, initialDc, approxLossDerivative, 0.001, 0.95, true, minSteps, maxSteps, func(dc float64) float64 {
		gammaLoss, gammaIntegral, gammaAnalytic, lastStepGammaConfidenceInterval, finalModel = gammaCalculationStep(itp, dc, parameters)
		gammaI = append(gammaI, gammaIntegral)
		return gammaLoss
	})
	gammaIntegral = utils.Mean(gammaI[max(len(gammaI)-minSteps, len(gammaI)/4):])
	gammaMargin = utils.StudentedMarginFromData(0.95, gammaI[max(len(gammaI)-minSteps, len(gammaI)/4):])
	gammaLoss = gammaAnalytic - gammaIntegral

	println("saved d_c and γ")
	return &GammaDataRow{
		OptimizationResult: OptimizationResult{ReducedFieldAtCathode: 2 * parameters.CathodeFallPotential / (dc * parameters.GasDensity * constants.Townsend),
			ReducedFieldAtSheathCenter: parameters.CathodeFallPotential / (dc * parameters.GasDensity * constants.Townsend),
			CathodeFallLength:          dc,
			GammaLoss:                  gammaLoss,
		},
		GammaAnalytic: gammaAnalytic,
		GammaIntegral: gammaIntegral,
		GammaMargin:   max(gammaMargin, lastStepGammaConfidenceInterval),
	}, finalModel
}
