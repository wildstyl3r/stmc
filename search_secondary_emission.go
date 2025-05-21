package main

import (
	"fmt"
	"math"

	"github.com/wildstyl3r/lxgata"
)

type LossType int

const (
	MSE LossType = iota
	Difference
)

var ionDriftVelocity = map[string](func(V, dc, N float64) float64){
	"Ar": func(V, dc, N float64) float64 { // Wald 1962
		E := V / dc
		return 4 * E / math.Cbrt(1.+math.Pow(0.007*E, 1.5))
	},
	"He": func(V, dc, N float64) float64 { // Raizer, Allen 1997
		E := V / (dc * N * Townsend)
		return 24. * E / math.Cbrt(1.+math.Pow(0.01*E, 1.5))
	},
}

func SnApproxIntegral(S0, lambda, dc, L float64) float64 {
	return -S0 * lambda * (math.Exp(-(L-dc)/lambda) - 1.)
}

func findXNumberDensityMaximum(de *DataExtractor) (xMaximum float64) {
	dcIndex := argmax(de.collisions[string(lxgata.IONIZATION)])
	dc := de.model.xStep * float64(dcIndex)
	if de.collisions[string(lxgata.IONIZATION)] != nil {
		S0 := de.collisions[string(lxgata.IONIZATION)][dcIndex] //* de.model.parameters.Pressure //* Torr / de.model.parameters.Pressure / de.cathodeFlux
		L := de.model.parameters.GapLength
		var lambda float64
		{
			SrInt := 0.
			for i := dcIndex; i < len(de.collisions[string(lxgata.IONIZATION)]); i++ {
				SrInt += de.collisions[string(lxgata.IONIZATION)][i]
			}
			SrInt *= de.model.xStep //* de.model.parameters.Pressure //* Torr / de.model.parameters.Pressure / de.cathodeFlux
			lambda = ternarySearchMax(func(lambda float64) float64 {
				analyticInt := SnApproxIntegral(S0, lambda, dc, L)
				return -(analyticInt - SrInt) * (analyticInt - SrInt)
			}, 0.001, 1000, 1.e-8)
		}
		//K := S0 * lambda * lambda / de.model.parameters.DAmbipolar //not important since we only need maximum
		if de.model.parameters.ParallelPlaneHollowCathode {
			xMaximum = L
		} else {
			xMaximum = ternarySearchMax(func(x float64) float64 {
				return glowDischargeNumberDensity(x, dc, L, lambda) //, D_amb float64)
			}, dc, L, 1.e-8)
		}
		return xMaximum
	} else {
		return math.NaN()
	}
}

func gammaIntegralF(de *DataExtractor) float64 {
	//calculate xMaximum
	xMaximum := findXNumberDensityMaximum(de)
	indexOfMaximum := min(int(xMaximum/de.model.xStep), len(de.collisions[string(lxgata.IONIZATION)]))
	// print("x_max/step: ", int(x_max/de.model.xStep), " len: ", len(de.collisions[string(lxgata.IONIZATION)]), "\n")
	sourceTermIntegral := de.collisions[string(lxgata.IONIZATION)][0]
	for i := 1; i < indexOfMaximum; i++ {
		sourceTermIntegral += de.collisions[string(lxgata.IONIZATION)][i]
	}
	sourceTermIntegral *= de.model.xStep * de.model.parameters.Pressure //* Torr / de.model.parameters.Pressure / de.cathodeFlux
	return 1 / sourceTermIntegral                                       // -> NaN -??

}

func gammaAnalyticF(dc, j, Vc, N float64, ionDriftVelocity func(float64, float64, float64) float64) float64 {
	return j*dc*dc/(2.*Vc*ionDriftVelocity(Vc, dc, N)*freeSpacePermittivityE0) - 1. //-> 0
}

func (dataExtractor *DataExtractor) gammaWithLoss(loss LossType) (lossValue float64, gi float64, ga float64) {
	ga = gammaAnalyticF(
		dataExtractor.model.parameters.CathodeFallLength,
		dataExtractor.model.parameters.CathodeCurrentDensity,
		dataExtractor.model.parameters.CathodeFallPotential,
		dataExtractor.model.parameters.gasDensity,
		ionDriftVelocity[dataExtractor.model.parameters.Species])
	gi = gammaIntegralF(dataExtractor)
	switch loss {
	case MSE:
		lossValue = (ga - gi) * (ga - gi)
	case Difference:
		lossValue = ga - gi
	}
	return
}

func estimateCathodeFallLengthLimits(parameters *ModelParameters) (from float64, to float64) {
	from, _ = binarySearch(func(dc float64) bool {
		return 0. < gammaAnalyticF(dc,
			parameters.CathodeCurrentDensity,
			parameters.CathodeFallPotential,
			parameters.gasDensity,
			ionDriftVelocity[parameters.Species])
	}, 1e-8, parameters.GapLength, parameters.CathodeFallLengthPrecision*0.01)

	_, to = binarySearch(func(dc float64) bool {
		return 1 < gammaAnalyticF(dc,
			parameters.CathodeCurrentDensity,
			parameters.CathodeFallPotential,
			parameters.gasDensity,
			ionDriftVelocity[parameters.Species]) // might be false everywhere in the gap, but as close as possible to true domain
	}, 0, parameters.GapLength, parameters.CathodeFallLengthPrecision*0.01)
	return from, to
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

func glowDischargeNumberDensity(x, d, L, lambda float64) float64 {
	return -(math.Exp(-(x-d)/lambda) - 1 - (math.Exp(-(L-d)/lambda)-1)*(x-d)/(L-d))
}

func gammaCalculationStep(itp *int, dc float64, parameters ModelParameters, lossType LossType) (loss, gammaIntegral, gammaAnalytic float64, dataExtractor *DataExtractor) {
	if itp != nil {
		fmt.Printf("step %d\n", *itp)
		*itp += 1
	} else {
		fmt.Println("preliminary step")
	}
	model := newModel(dc, parameters)
	model.run()
	dataExtractor = newDataExtractor(&model)
	loss, gammaIntegral, gammaAnalytic = dataExtractor.gammaWithLoss(lossType)
	if parameters._verbose {
		fmt.Printf("d_c: %v\nsecondary emission coefficient\n\t integral: %6f\n\t analytic:%6f\n", dc, gammaIntegral, gammaAnalytic)
	}
	return
}
