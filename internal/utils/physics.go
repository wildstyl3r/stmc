package utils

import (
	"math"

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

var SimplifiedIonDriftVelocity = map[string](func(V, dc, N float64) float64){
	"Ar": func(V, dc, N float64) float64 { // Wald 1962
		E := V / (dc * N * constants.Townsend)
		return 47.809144373375744 * math.Sqrt(E)
	},
	"He": func(V, dc, N float64) float64 { // Raizer, Allen 1997
		E := V / (dc * N * constants.Townsend)
		return 240 * math.Sqrt(E)
	},
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
