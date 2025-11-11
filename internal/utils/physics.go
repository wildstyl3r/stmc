package utils

import (
	"math"

	"github.com/wildstyl3r/lxgata"
	"github.com/wildstyl3r/stmc/internal/constants"
)

var IonDriftVelocity = map[string](func(V, dc, N float64) float64){
	"Ar": func(V, dc, N float64) float64 { // Wald 1962
		E := V / (dc * N * constants.Townsend)
		return 4 * E / math.Cbrt(1.+math.Pow(0.007*E, 1.5))
	},
	"He": func(V, dc, N float64) float64 { // Raizer, Allen 1997
		E := V / (dc * N * constants.Townsend)
		return 24. * E / math.Cbrt(1.+math.Pow(0.01*E, 1.5))
	},
}

var AtomicNumbers = map[string]lxgata.AtomicNumber{
	"He": 2,
	"Ar": 18,
}

var OpalIonizationShapeParameters = map[string]float64{
	"He": 15.8,
	"Ne": 24.2,
	"H2": 8.3,
	"N2": 13.0,
	"O2": 17.4,
	"CO": 14.2,
	"NO": 13.6,

	"Ar": 10, // Ar, Kr and Xe are marked as unsure in the original paper (doi: 10.1063/1.1676707)
	"Kr": 9.6,
	"Xe": 8.7,
}

var WeakFieldConstantIonMobilityPerTorr = map[string]float64{ // mu/p[Torr]
	"Ar": 0.0833,
	"He": 0.83,
}

var SimplifiedIonDriftVelocityCoefficient = map[string]float64{ // v_{id}(E/n) ~ K*(E/n)^0.5
	"Ar": 47.809144373375744,
	"He": 240,
}

func EV2J(val float64) float64 {
	return val * constants.ElementaryCharge
}

func J2eV(val float64) float64 {
	return val / constants.ElementaryCharge
}

func EV2electronVelocity(energy float64) (v float64) {
	v = math.Sqrt(2. * EV2J(energy) / constants.ElectornMass)
	return
}
