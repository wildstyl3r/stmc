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
