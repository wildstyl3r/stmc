package main

import (
	"math"
	"math/rand"
)

const k float64 = 1.380649e-23
const electronCharge = 1.602176634e-19 // C
const me float64 = 9.1093837139e-31    // [kg]
const Torr float64 = 101325. / 760.    // [Pa]

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

type AveragingElement struct {
	sum    float64
	number int
}

func (a *AveragingElement) avg() float64 {
	if a.number == 0 {
		return 0
	}
	return a.sum / float64(a.number)
}

func (a *AveragingElement) add(x float64) {
	a.sum += x
	a.number += 1
}

type StateIncrement struct {
	xIndex, eIndex, muIndex int
}

type ProbIncrement struct {
	x        int
	collType string
	prob     float64
}

type CollisionEvent struct {
	x          int
	energyLoss float64
	collType   string
}
