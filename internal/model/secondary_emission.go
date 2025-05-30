package model

import (
	"fmt"
	"math"

	"github.com/wildstyl3r/lxgata"
	"github.com/wildstyl3r/stmc/internal/config"
	"github.com/wildstyl3r/stmc/internal/constants"
	"github.com/wildstyl3r/stmc/internal/utils"
)

var IonDriftVelocity = map[string](func(V, dc, N float64) float64){
	"Ar": func(V, dc, N float64) float64 { // Wald 1962
		E := V / dc
		return 4 * E / math.Cbrt(1.+math.Pow(0.007*E, 1.5))
	},
	"He": func(V, dc, N float64) float64 { // Raizer, Allen 1997
		E := V / (dc * N * constants.Townsend)
		return 24. * E / math.Cbrt(1.+math.Pow(0.01*E, 1.5))
	},
}

func GlowDischargeDensity(de *DataExtractor) (density []float64) {
	density = make([]float64, len(de.collisions[lxgata.IONIZATION]))

	dcIndex := min(int(de.model.Parameters.CathodeFallLength/de.model.XStep), de.model.NumCells-2)
	exactSnAtDc := de.collisions[lxgata.IONIZATION][dcIndex] + (de.model.Parameters.CathodeFallLength-float64(dcIndex)*de.model.XStep)/de.model.XStep*(de.collisions[lxgata.IONIZATION][dcIndex+1]-de.collisions[lxgata.IONIZATION][dcIndex])
	C1 := 0.5 * (exactSnAtDc + de.collisions[lxgata.IONIZATION][dcIndex]) * de.model.XStep
	{
		accum := utils.TableIntegrate(de.collisions[lxgata.IONIZATION][:dcIndex], nil, de.model.XStep)
		for i := dcIndex; i < len(de.collisions[lxgata.IONIZATION]); i++ {
			C1 += accum
			accum += de.collisions[lxgata.IONIZATION][i] * de.model.XStep
		}
		C1 *= -de.model.XStep
		C1 /= (de.model.Parameters.GapLength - de.model.Parameters.CathodeFallLength)
	}
	C2 := 0.
	{
		accum := 0.
		for i := range de.collisions[lxgata.IONIZATION] {
			C2 -= accum
			accum += de.collisions[lxgata.IONIZATION][i] * de.model.XStep //TableIntegrate(de.collisions[lxgata.IONIZATION][:i], nil, de.model.XStep)
		}
		C2 *= de.model.XStep
		C2 -= C1 * de.model.Parameters.GapLength
	}
	{
		sum := 0.
		accum := 0.
		for j := range de.collisions[lxgata.IONIZATION] {
			sum += accum
			accum += de.collisions[lxgata.IONIZATION][j] * de.model.XStep
			density[j] = -(sum*de.model.XStep + C1*float64(j)*de.model.XStep + C2)
		}
	}
	return
}

func gammaIntegralF(de *DataExtractor) float64 {
	//calculate xMaximum
	// xMaximum := findXNumberDensityMaximum(de)
	indexOfMaximum := utils.Argmax(GlowDischargeDensity(de))
	// print("x_max/step: ", int(x_max/de.model.XStep), " len: ", len(de.collisions[string(lxgata.IONIZATION)]), "\n")
	sourceTermIntegral := utils.TableIntegrate(de.collisions[lxgata.IONIZATION][0:indexOfMaximum], nil, de.model.XStep) //0. //de.collisions[lxgata.IONIZATION][0]
	sourceTermIntegral *= de.model.Parameters.Pressure                                                                  //* Torr / de.model.Parameters.Pressure / de.cathodeFlux
	return 1 / sourceTermIntegral                                                                                       // -> NaN -??

}

func gammaAnalyticF(dc, j, Vc, N float64, ionDriftVelocity func(float64, float64, float64) float64) float64 {
	return j*dc*dc/(2.*Vc*ionDriftVelocity(Vc, dc, N)*constants.FreeSpacePermittivityE0) - 1. //-> 0
}

func (dataExtractor *DataExtractor) gammaWithLoss(loss utils.LossType) (lossValue float64, gi float64, ga float64) {
	ga = gammaAnalyticF(
		dataExtractor.model.Parameters.CathodeFallLength,
		dataExtractor.model.Parameters.CathodeCurrentDensity,
		dataExtractor.model.Parameters.CathodeFallPotential,
		dataExtractor.model.Parameters.GasDensity,
		IonDriftVelocity[dataExtractor.model.Parameters.Species])
	gi = gammaIntegralF(dataExtractor)
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
			IonDriftVelocity[parameters.Species])
	}, 1e-8, parameters.GapLength, parameters.CathodeFallLengthPrecision*0.01)

	_, to = utils.BinarySearch(func(dc float64) bool {
		return 1 < gammaAnalyticF(dc,
			parameters.CathodeCurrentDensity,
			parameters.CathodeFallPotential,
			parameters.GasDensity,
			IonDriftVelocity[parameters.Species]) // might be false everywhere in the gap, but as close as possible to true domain
	}, 0, parameters.GapLength, parameters.CathodeFallLengthPrecision*0.01)
	return from, to
}

func GetApproximateDcForGamma(gamma, minDc, maxDc float64, parameters config.ModelParameters) float64 {
	initialDcL, initialDcR := utils.BinarySearch(func(dc float64) bool {
		return gamma < gammaAnalyticF(
			dc,
			parameters.CathodeCurrentDensity,
			parameters.CathodeFallPotential,
			parameters.GasDensity,
			IonDriftVelocity[parameters.Species])
	}, minDc, maxDc, parameters.CathodeFallLengthPrecision*0.1)
	return 0.5 * (initialDcL + initialDcR)
}

// func delta_plus(L, d, lambda, big_lambda float64) float64 {
// 	return (1 + math.Exp(-(L-2*d)/lambda)) / (1 + math.Exp(-(L-2*d)/big_lambda))
// }

// func pphc_n(x, d, L, lambda, big_lambda float64) float64 { //, D_amb
// 	// lambda2 := lambda * lambda
// 	// return lambda2 / (D_amb * (1 - lambda2/(big_lambda*big_lambda))) *
// 	// (S0 / (1 - math.Exp(-(L-2*d)/lambda))) *
// 	return delta_plus(L, d, lambda, big_lambda)*(math.Exp(-(x-d)/big_lambda)+math.Exp(-(L-x-d)/big_lambda)) - math.Exp(-(x-d)/lambda) - math.Exp(-(L-x-d)/lambda)
// }

func GammaCalculationStep(itp *int, dc float64, parameters config.ModelParameters, lossType utils.LossType) (loss, gammaIntegral, gammaAnalytic float64, dataExtractor *DataExtractor) {
	if itp != nil {
		fmt.Printf("step %d\n", *itp)
		*itp += 1
	} else {
		fmt.Println("preliminary step")
	}
	model := NewModel(dc, parameters)
	model.Run()
	dataExtractor = NewDataExtractor(&model)
	loss, gammaIntegral, gammaAnalytic = dataExtractor.gammaWithLoss(lossType)
	if parameters.Verbose() {
		fmt.Printf("d_c: %v\nsecondary emission coefficient\n\t integral: %6f\n\t analytic:%6f\n", dc, gammaIntegral, gammaAnalytic)
	}
	return
}
