package utils

import (
	"math"

	"golang.org/x/exp/constraints"
	"gonum.org/v1/gonum/stat/distuv"
)

type Number interface {
	constraints.Float | constraints.Integer
}

func Mean[T Number](s []T) (mean float64) {
	floated := make([]float64, len(s))
	for i := range s {
		floated[i] = float64(s[i])
	}
	mean = SumFloat64Slice(floated) / float64(len(s))
	return
}

func MeanAndVariance[T Number](s []T, unbiased bool) (mean, variance float64) {
	mean = Mean(s)
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

func NormalMargin(quantileCoefficient, variance float64, sampleSize int) float64 {
	return quantileCoefficient * math.Sqrt(variance/float64(sampleSize))
}

func UnbiasedThirdCentralMomentEstimation(sample []float64) float64 {
	n := float64(len(sample))
	sumTerms := make([]float64, len(sample))
	mean := Mean(sample)
	for i := range sample {
		sumTerms[i] = math.Pow(sample[i]-mean, 3)
	}
	return n * SumFloat64Slice(sumTerms) / ((n - 1) * (n - 2))
}

func RawMomentsFrom1UpTo(order int, sample []float64) (raw []float64) {
	raw = make([]float64, order)
	n := float64(len(sample))
	for r := 1; r <= order; r++ {
		sumTerms := make([]float64, len(sample))
		for i := range sample {
			sumTerms[i] = math.Pow(sample[i], float64(r))
		}
		raw[r-1] = SumFloat64Slice(sumTerms) / n
	}
	return raw
}

func SecondOrderTaylorInverseExpectationEstimation(rawMoments []float64) float64 {
	return rawMoments[1] / math.Pow(rawMoments[0], 3)
}

func SecondOrderTaylorInverseVarianceEstimation(rawMoments []float64) float64 {
	A, B, C, D := rawMoments[0], rawMoments[1], rawMoments[2], rawMoments[3]
	sumTerms := []float64{
		3.,
		-3 * B / (A * A),
		(D - B*B) / math.Pow(A, 6),
		-6 * C / math.Pow(A, 5),
		6 * B / math.Pow(A, 4),
	}
	return SumFloat64Slice(sumTerms)
}

func StudentedMargin(confidenceP float64, data []float64) float64 {
	dist := distuv.StudentsT{Nu: float64(len(data) - 1), Sigma: 1}
	quantileFactor := dist.Quantile(1 - (1-confidenceP)/2) // make it two-sided
	return quantileFactor * math.Sqrt(Variance(data, true)/float64(len(data)))
}
