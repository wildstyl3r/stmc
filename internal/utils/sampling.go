package utils

import (
	"math"
	"math/rand/v2"
)

func R() float64 {
	return rand.ExpFloat64()
}

func UniformOnDisk(r float64) (a, b float64) {
	a, b = 2.*rand.Float64()-1., 2.*rand.Float64()-1.
	for a*a+b*b > 1. {
		a, b = 2.*rand.Float64()-1., 2.*rand.Float64()-1.
	}
	a *= r
	b *= r
	return
}

func BoxMuller2Normals() (n1, n2 float64) {
	x, y := UniformOnDisk(1)
	s := x*x + y*y
	f := math.Sqrt(-2 * math.Log(s) / s)
	return x * f, y * f
}
