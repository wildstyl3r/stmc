package extensions

import (
	"cmp"
	"fmt"
	"math"
	"os"

	"github.com/aclements/go-moremath/fit"
	"github.com/schollz/progressbar/v3"
	"github.com/wildstyl3r/stmc/internal/config"
	"github.com/wildstyl3r/stmc/internal/model"
	"github.com/wildstyl3r/stmc/internal/utils"
)

func sourceIntegralMonteCarloF(parameters config.ModelParameters, requirePrecision bool) (sumMean, sumMargin float64, m *model.Model) {
	m = model.NewModel(parameters)
	LoadExtensions(m.DataHub)
	m.Run(func(m *model.Model) int {
		if m.TotalElectronsEmittedOnCathode == 0 {
			return m.Parameters.NElectrons
		} else {
			sumMean, sumMargin = m.GetMetrics(SourceIntegralPerCathodeElectronFluxKey).(float64), m.GetMetrics(SourceIntegralPerCathodeElectronFluxMarginKey).(float64)
			if parameters.SourceIntegralRelativeMargin < sumMargin/sumMean && parameters.SourceIntegralRelativeMargin != 0 && requirePrecision {
				return m.Parameters.AddByNElectrons
			} else {
				return 0
			}
		}
	})
	if math.IsNaN(sumMean) {
		fmt.Println("integral is NaN")
	}
	return sumMean, sumMargin, m
}

// not applicable for UniformField
func sourceIntegralAnalyticPhelpsF(parameters config.ModelParameters) float64 {
	return 1. / gammaAnalyticPhelpsF(parameters) //-> 0
}

func gammaDonkoF(parameters config.ModelParameters) float64 {
	return 0.01 * math.Pow(parameters.ReducedFieldAtCathode()/1000, 0.6)
}

func sourceIntegralAnalyticF(parameters config.ModelParameters) float64 {
	return 1. / gammaAnalyticF(parameters)
}
func gammaAnalyticF(parameters config.ModelParameters) float64 {
	if parameters.CalculateSecondaryEmissionCoefficient {
		return gammaAnalyticPhelpsF(parameters)
	} else if parameters.CalculateCurrentDensity || parameters.CalculateVoltage {
		if parameters.UseDonkosSEC {
			return gammaDonkoF(parameters)
		} else if parameters.SecondaryEmissionCoefficient != 0 {
			return parameters.SecondaryEmissionCoefficient
		} else {
			panic("unable to determine analytic source integral")
		}
	}
	return 0
}

// func getApproximateCathodeFallLengthForSourceIntegral(sourceIntegral, minDc, maxDc float64, parameters config.ModelParameters) float64 {
// 	initialDcL, initialDcR := utils.BinarySearch(func(dc float64) bool {
// 		parameters.CathodeFallLength = dc
// 		si := sourceIntegralAnalyticPhelpsF(parameters)
// 		return sourceIntegral > si
// 	}, minDc, maxDc, parameters.CathodeFallLengthPrecision*0.001)
// 	return 0.5 * (initialDcL + initialDcR)
// }

func sourceIntegralCalculationStep(requirePrecision bool, parameters *config.ModelParameters, analyticSourceIntegral func(config.ModelParameters) float64) (optResult OptimizationResult, sourceIntegralVariance float64, stepModel *model.Model) {
	sourceIntegralMonteCarlo, sourceIntegralVariance, stepModel := sourceIntegralMonteCarloF(*parameters, requirePrecision)
	sourceIntegralAnalytic := analyticSourceIntegral(*parameters) //sourceIntegralAnalyticF(stepModel.Parameters.CathodeFallLength, stepModel.Parameters.CathodeCurrentDensity, stepModel.Parameters.CathodeFallPotential, stepModel.Parameters.GasDensity, utils.IonDriftVelocity[stepModel.Parameters.Species])
	difference := sourceIntegralAnalytic - sourceIntegralMonteCarlo
	return OptimizationResult{
			ModelName:                               parameters.PrototypeName(),
			ReducedFieldAtCathode:                   stepModel.ReducedFieldAtCathode(),
			ReducedFieldAtSheathCenter:              stepModel.ReducedFieldMidSheath(),
			CathodeFallLength:                       stepModel.Parameters.CathodeFallLength,
			PressureCathodeFallLength:               stepModel.Parameters.CathodeFallLength * stepModel.Parameters.Pressure,
			SourceIntegralDifference:                difference,
			SourceIntegralAnalytic:                  sourceIntegralAnalytic,
			SourceIntegralMonteCarlo:                sourceIntegralMonteCarlo,
			SourceIntegralMargin:                    utils.EstimateMargin(0.95, sourceIntegralVariance, float64(stepModel.TotalElectronsEmittedOnCathode)),
			GammaAnalytic:                           1. / sourceIntegralAnalytic,
			GammaMonteCarlo:                         1. / sourceIntegralMonteCarlo,
			MeanElectronEnergyAtAnode:               stepModel.MeanElectronEnergyAtAnode,
			MeanFreePathAnode:                       1. / (stepModel.Parameters.CrossSectionsData().TotalCrossSectionAt(stepModel.MeanElectronEnergyAtAnode) * stepModel.Parameters.GasDensity),
			GlobalMeanFreePath:                      stepModel.GlobalMeanFreePath.Mean(),
			Voltage:                                 stepModel.Parameters.CathodeFallPotential,
			Pressure:                                stepModel.Parameters.Pressure,
			CathodeCurrent:                          stepModel.Parameters.CathodeCurrentDensity * (math.Pi * parameters.CathodeRadius * parameters.CathodeRadius),
			CathodeCurrentDensity:                   stepModel.Parameters.CathodeCurrentDensity,
			CathodeCurrentDensityPerPressureSquared: stepModel.Parameters.CathodeCurrentDensity / (parameters.Pressure * parameters.Pressure),
		},
		sourceIntegralVariance,
		stepModel
}

type SourceIntegralDataRow struct {
	OptimizationResult
	GammaMonteCarloMargin float64 `csv:"$\\gamma$ Monte Carlo margin"`
}

func SourceIntegralCalculation(parameters config.ModelParameters, outputDir, modelName string) (dataRow SourceIntegralDataRow, altDataRow *SourceIntegralDataRow, finalmodel, altModel *model.Model) {
	_, maxDc := EstimateCathodeFallLengthLimits(parameters)
	minDc := parameters.CathodeFallLengthPrecision
	numberOfSteps := int((maxDc - minDc) / parameters.CathodeFallLengthPrecision)
	return GeneralizedCalculation(parameters, numberOfSteps, minDc, parameters.CathodeFallLengthPrecision,
		func(dc float64, mp *config.ModelParameters) {
			mp.CathodeFallLength = dc
		},
		func(optResult OptimizationResult, variance float64, m *model.Model) SourceIntegralDataRow {
			return SourceIntegralDataRow{
				OptimizationResult:    optResult,
				GammaMonteCarloMargin: (math.Abs(optResult.SourceIntegralDifference) + optResult.SourceIntegralMargin) / (optResult.SourceIntegralMonteCarlo * optResult.SourceIntegralMonteCarlo),
			}
		}, func(p config.ModelParameters, simc, sia float64) bool {
			return (sia > 0 && simc > 3*sia)
		}, true, outputDir, modelName)
}

func GeneralizedCalculation[K cmp.Ordered, T utils.Indexable[K]](parameters config.ModelParameters,
	numberOfSteps int, offset, stepSize float64,
	setUp func(arg float64, mp *config.ModelParameters),
	extractor func(optResult OptimizationResult, variance float64, m *model.Model) T,
	stopCondition func(p config.ModelParameters, simc, sia float64) bool,
	preferGreaterRoot bool,
	outputDir, modelName string,
) (dataRow T, altDataRow *T, finalModel, altModel *model.Model) {
	//imprecise monte-carlo data gathering
	//regression fitting

	//find root on source difference

	//return precise monte-carlo calculation on root
	stepArg, stepSourceIntegralMC, stepVariance := make([]float64, 0, numberOfSteps), make([]float64, 0, numberOfSteps), make([]float64, 0, numberOfSteps)
	minArg, maxArg := offset, offset+float64(numberOfSteps)*stepSize
	{
		var debugData []T = make([]T, 0, numberOfSteps)
		progress := progressbar.Default(int64(numberOfSteps))
		parameters.SupressSpinner = true
		for i := range numberOfSteps {

			// fmt.Printf("step %d of %d\n", i+1, numberOfSteps)
			arg := offset + float64(i)*stepSize
			setUp(arg, &parameters)
			optResult, variance, m := sourceIntegralCalculationStep(false, &parameters, sourceIntegralAnalyticF)
			progress.Describe(fmt.Sprintf("[ analytic: %6f ; Monte Carlo: %6f +/- %.2f%% ]", optResult.SourceIntegralAnalytic, optResult.SourceIntegralMonteCarlo, optResult.SourceIntegralMargin/optResult.SourceIntegralMonteCarlo*100))
			stepArg, stepSourceIntegralMC, stepVariance = append(stepArg, arg), append(stepSourceIntegralMC, optResult.SourceIntegralMonteCarlo), append(stepVariance, variance)
			if parameters.DebugOutput {
				row := extractor(optResult, variance, m)
				debugData = append(debugData, row)
			}
			if len(stepArg) > 20 && stopCondition(parameters, optResult.SourceIntegralMonteCarlo, optResult.SourceIntegralAnalytic) {
				maxArg = offset + float64(i)*stepSize
				// progress.Add(numberOfSteps - i)
				fmt.Println("\nearly stop condition is met")
				break
			}
			progress.Add(1)
		}
		parameters.SupressSpinner = false
		if parameters.DebugOutput {
			err := utils.WriteAsCSV(config.MakeHeader(debugData[0], parameters.OutputUnits()), debugData, outputDir+"/"+modelName, "source_integral_bruteforce", false)

			if err != nil {
				println("unable to save dc and secondary emission coefficient: ", err.Error())
				os.Exit(1)
			}
		}
	}
	regression := fit.PolynomialRegression(stepArg, stepSourceIntegralMC, utils.Apply(func(x float64) float64 { return 1. / x }, stepVariance), 2)
	fmt.Printf("Regression coefficients are %v\n", regression.Coefficients)
	var sourceDifferenceRoot float64
	altSourceDifferenceRoot := math.Inf(-1)
	{

		_, minArg = utils.BinarySearch(func(arg float64) bool {
			setUp(arg, &parameters)
			return sourceIntegralAnalyticF(parameters) > 0
		}, minArg, maxArg, stepSize*0.00001)

		equation := func(arg float64) float64 {
			setUp(arg, &parameters)
			analytic := sourceIntegralAnalyticF(parameters)
			approx := regression.F(arg)
			return analytic - approx
		}
		if math.Signbit((equation(minArg+stepSize)-equation(minArg))/stepSize) == math.Signbit((equation(maxArg)-equation(maxArg-stepSize))/stepSize) {
			baseline := equation(minArg) < 0
			left, right := utils.BinarySearch(func(arg float64) bool {
				return equation(arg) < 0 != baseline
			}, minArg, maxArg, stepSize*0.000001)
			sourceDifferenceRoot = 0.5 * (left + right)
			fmt.Printf("The only root is at %f, eq at root: %f\n", sourceDifferenceRoot, equation(sourceDifferenceRoot))
		} else {
			var maxDiff float64
			if (equation(minArg+stepSize)-equation(minArg))/stepSize > 0 {
				maxDiff = utils.TernarySearchMax(func(arg float64) float64 {
					return equation(arg)
				}, minArg, maxArg, stepSize*0.000001)

			} else {
				maxDiff = utils.TernarySearchMax(func(arg float64) float64 {
					return -equation(arg)
				}, minArg, maxArg, stepSize*0.000001)
			}
			baseline := equation(minArg) < 0
			firstRootL, firstRootR := utils.BinarySearch(func(arg float64) bool {
				return equation(arg) < 0 != baseline
			}, minArg, maxDiff, stepSize*0.000001)
			firstRoot := 0.5 * (firstRootL + firstRootR)

			baseline = equation(maxArg) < 0
			secondRootL, secondRootR := utils.BinarySearch(func(arg float64) bool {
				return equation(arg) < 0 == baseline
			}, maxDiff, maxArg, stepSize*0.000001)
			secondRoot := 0.5 * (secondRootL + secondRootR)

			if math.Abs(firstRoot-secondRoot) < stepSize {
				sourceDifferenceRoot = 0.5 * (firstRoot + secondRoot)
				fmt.Printf("No exact root. The closest point is %f with difference %f\n", sourceDifferenceRoot, equation(sourceDifferenceRoot))
			} else if math.Abs(equation(firstRoot)) > 1e-3 || math.Abs(equation(secondRoot)) > 1e-3 {
				if math.Abs(equation(firstRoot)) < math.Abs(equation(secondRoot)) {
					sourceDifferenceRoot = firstRoot
				} else {
					sourceDifferenceRoot = secondRoot
				}
				fmt.Printf("At least one of the roots is outside the boundaries, the best accepted: %f, eq at accepted: %f\n", sourceDifferenceRoot, equation(sourceDifferenceRoot))
			} else {
				fmt.Printf("First root %f, eq at first root is %f, second root %f, eq at second root is %f\n", firstRoot, equation(firstRoot), secondRoot, equation(secondRoot))
				if preferGreaterRoot {
					sourceDifferenceRoot = secondRoot
					altSourceDifferenceRoot = firstRoot
				} else {
					sourceDifferenceRoot = firstRoot
					altSourceDifferenceRoot = secondRoot
				}
			}
		}
	}
	setUp(sourceDifferenceRoot, &parameters)
	optResult, variance, m := sourceIntegralCalculationStep(true, &parameters, sourceIntegralAnalyticF)
	fmt.Printf("[ analytic: %6f ; Monte Carlo: %6f +/- %.2f%% ]\n", optResult.SourceIntegralAnalytic, optResult.SourceIntegralMonteCarlo, optResult.SourceIntegralMargin/optResult.SourceIntegralMonteCarlo*100)
	optResult.CathodeFallLengthMargin = math.Abs(optResult.SourceIntegralMargin / utils.PolynomialDerivative(optResult.CathodeFallLength, regression.Coefficients))
	optResult.PressureCathodeFallLengthMargin = optResult.CathodeFallLengthMargin * optResult.Pressure
	dataRow = extractor(optResult, variance, m)
	if !math.IsInf(altSourceDifferenceRoot, -1) {
		setUp(altSourceDifferenceRoot, &parameters)
		altResult, altVariance, altM := sourceIntegralCalculationStep(true, &parameters, sourceIntegralAnalyticF)
		fmt.Printf("[ analytic: %6f ; Monte Carlo: %6f +/- %.2f%% ]\n", altResult.SourceIntegralAnalytic, altResult.SourceIntegralMonteCarlo, altResult.SourceIntegralMargin/altResult.SourceIntegralMonteCarlo*100)
		altResult.CathodeFallLengthMargin = math.Abs(altResult.SourceIntegralMargin / utils.PolynomialDerivative(altResult.CathodeFallLength, regression.Coefficients))
		altResult.PressureCathodeFallLengthMargin = altResult.CathodeFallLengthMargin * altResult.Pressure
		altDataRow := extractor(altResult, altVariance, altM)
		return dataRow, &altDataRow, m, altM
	} else {
		return dataRow, nil, m, nil
	}
}
