package main

import (
	"math"

	"github.com/wildstyl3r/lxgata"
)

func exp_integral(S0, lambda, dc, L float64) float64 {
	return -S0 * lambda * (math.Exp(-(L-dc)/lambda) - 1.)
}

func find_x_max(de *DataExtractor) (x_max float64) {
	dc_i := argmax(de.collisions[string(lxgata.IONIZATION)])
	dc := de.model.xStep * float64(dc_i)
	S0 := de.collisions[string(lxgata.IONIZATION)][dc_i] //* Torr / de.model.parameters.Pressure / de.cathodeFlux
	L := de.model.parameters.GapLength
	var lambda float64
	{

		SrInt := 0.
		for i := dc_i; i < len(de.collisions[string(lxgata.IONIZATION)]); i++ {
			SrInt += de.collisions[string(lxgata.IONIZATION)][i]
		}
		SrInt *= de.model.xStep //* Torr / de.model.parameters.Pressure / de.cathodeFlux
		lambda = ternarySearchMax(func(lambda float64) float64 { return exp_integral(S0, lambda, dc, L) }, 0.001, 1000, 1.e-8)
	}
	//K := S0 * lambda * lambda / de.model.parameters.DAmbipolar //not important since we only need maximum
	if de.model.parameters.ParallelPlaneHollowCathode {
		x_max = L
	} else {
		x_max = ternarySearchMax(func(x float64) float64 {
			return glow_discharge_n(x, dc, L, lambda) //, D_amb float64)
		}, dc, L, 1.e-8)
	}
	return x_max
}

func gamma_integral(de *DataExtractor) float64 {
	//calculate x_max
	x_max := find_x_max(de)
	i_max := min(int(x_max/de.model.xStep), len(de.collisions[string(lxgata.IONIZATION)]))
	S_integral := 0.
	for i := 1; i < i_max; i++ {
		S_integral += de.collisions[string(lxgata.IONIZATION)][i]
	}
	S_integral *= de.model.xStep //* Torr / de.model.parameters.Pressure / de.cathodeFlux
	return 1 / S_integral        // -> NaN -??

}

func gamma_analytic(de *DataExtractor) float64 {
	dc := de.model.parameters.CathodeFallLength
	j := de.model.parameters.CathodeCurrentDensity
	Vc := de.model.parameters.CathodeFallPotential
	v_i_d := v_i_d_Helium_Raizer(Vc, dc, de.model.gasDensity)
	return j * dc * dc / (2. * Vc * v_i_d * e0) //-> 0
}

func v_i_d_Helium_Raizer(V, dc, N float64) float64 {
	E := V / (dc * N * Townsend)
	return 24. * E / math.Cbrt(1.+math.Pow(0.01*E, 1.5))
}

func (dataExtractor *DataExtractor) gamma_loss() (float64, float64, float64) {
	gi := gamma_integral(dataExtractor)
	ga := gamma_analytic(dataExtractor)
	return (gi - ga) * (gi - ga), gi, ga
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

func glow_discharge_n(x, d, L, lambda float64) float64 {
	return -(math.Exp(-(x-d)/lambda) - 1 - (math.Exp(-(L-d)/lambda)-1)*(x-d)/(L-d))
}
