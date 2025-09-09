package utils

import (
	"math"
	"math/rand"
)

func R() float64 {
	return -math.Log(rand.Float64())
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
