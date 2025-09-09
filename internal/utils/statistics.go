package utils

import "golang.org/x/exp/constraints"

type Number interface {
	constraints.Float | constraints.Integer
}

func Average[T Number](s []T) (mean float64) {
	for i := range s {
		mean += float64(s[i])
	}
	mean /= float64(len(s))
	return
}

func MeanAndVariance[T Number](s []T, unbiased bool) (mean, variance float64) {
	mean = Average(s)
	for i := range s {
		variance += (float64(s[i]) - mean) * (float64(s[i]) - mean)
	}
	if unbiased {
		variance /= float64(len(s) - 1)
	} else {
		variance /= float64(len(s))
	}

	return
}

func Variance[T Number](s []T, unbiased bool) float64 {
	_, v := MeanAndVariance(s, unbiased)
	return v
}
