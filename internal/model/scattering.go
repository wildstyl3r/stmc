package model

import (
	"math"
	"math/rand"
)

type ScatteringFunction func(energy, energyAfter, energyLoss float64) (cos float64)

func isotropic(energy, energyAfter, energyLoss float64) float64 {
	return 1. - 2.*rand.Float64()
}

func surendra(energy, energyAfter, energyLoss float64) float64 {
	return (2. + energy - 2.*math.Pow(1.+energy, rand.Float64())) / energy
}

func bornDipole(energy, energyAfter, energyLoss float64) float64 {
	energyRatioSquare := energyLoss * energyLoss / math.Pow(math.Sqrt(energyAfter)+math.Sqrt(energy), 4)
	return 1. + 2.*energyRatioSquare/(1.-energyRatioSquare)*(1.-math.Pow(energyRatioSquare, -rand.Float64()))

}
