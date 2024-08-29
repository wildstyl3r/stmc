package main

import (
	"math"
	"math/rand"
)

func ternarySearchMaxF(f func(float64) float64, left, right, eps float64) float64 {
	for right-left > eps {
		a := (left*2. + right) / 3.
		b := (left + right*2.) / 3.
		if f(a) > f(b) {
			right = b
		} else {
			left = a
		}
	}
	return f((left + right) / 2.)
}

func R() float64 {
	return -math.Log(rand.Float64())
}

func cm2m(val float64) float64 {
	return val / 100.
}
func m2cm(val float64) float64 {
	return val * 100.
}

func Da2kg(val float64) float64 {
	return val * 1.66053906892e-27
}

func eV2J(val float64) float64 {
	return val * 1.602176634e-19
}

func J2eV(val float64) float64 {
	return val / 1.602176634e-19
}
