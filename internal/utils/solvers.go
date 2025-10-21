package utils

import "math"

func TernarySearchMax(f func(float64) float64, left, right, eps float64) float64 {
	for right-left > eps {
		a := math.FMA(left, 2., right) / 3.
		b := math.FMA(right, 2., left) / 3.
		if f(a) > f(b) {
			right = b
		} else {
			left = a
		}
	}
	return (left + right) * 0.5
}

// return the point of the condition support that is not farther than eps from the support boundary
// invariant: at *right* condition must be TRUE
func BinarySearch(condition func(float64) bool, falseDom, trueDom, eps float64) (float64, float64) {
	for math.Abs(trueDom-falseDom) > eps {
		c := (falseDom + trueDom) * 0.5
		if condition(c) {
			trueDom = c
		} else {
			falseDom = c
		}
	}
	return falseDom, trueDom
}

func TernarySearchMaxF(f func(float64) float64, left, right, eps float64) float64 {
	return f(TernarySearchMax(f, left, right, eps))
}

func IntersectContinuousWithQuasiLinear(lowerBound, upperBound, argEps, fEps float64,
	continuous func(float64) float64, quasiLinear func(float64) float64) (arg, val float64) {
	lowerQL, upperQL := quasiLinear(lowerBound), quasiLinear(upperBound)
	delta := math.Inf(1)
	for math.Abs(delta) > fEps && math.Abs(upperBound-lowerBound) > argEps {
		k := (upperQL - lowerQL) / (upperBound - lowerBound)
		b := lowerQL - k*lowerBound
		crossingEstimateF, crossingEstimateT := BinarySearch(func(x float64) bool {
			return continuous(x)-(k*x+b) < 0
		}, lowerBound, upperBound, argEps)
		crossingEstimate := 0.5 * (crossingEstimateF + crossingEstimateT)
		evalCrossingQL := quasiLinear(crossingEstimate)
		evalCrossingDC := continuous(crossingEstimate)
		delta = evalCrossingDC - evalCrossingQL
		if evalCrossingDC-evalCrossingQL > 0 {
			lowerBound = crossingEstimate
			lowerQL = evalCrossingQL
		} else {
			upperBound = crossingEstimate
			upperQL = evalCrossingQL
		}
	}
	return 0.5 * (lowerBound + upperBound), continuous(0.5 * (lowerBound + upperBound))
}
