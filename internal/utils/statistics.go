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

func StudentedMarginFromData(confidenceP float64, data []float64) float64 {
	dist := distuv.StudentsT{Nu: float64(len(data) - 1), Sigma: 1}
	quantileFactor := dist.Quantile(1 - (1-confidenceP)/2) // make it two-sided
	return quantileFactor * math.Sqrt(Variance(data, true)/float64(len(data)))
}
func StudentedMargin(confidenceP, variance float64, sampleSize int) float64 {
	dist := distuv.StudentsT{Nu: float64(sampleSize - 1), Sigma: 1}
	quantileFactor := dist.Quantile(1 - (1-confidenceP)/2) // make it two-sided
	return quantileFactor * math.Sqrt(variance/float64(sampleSize))
}

type SampleCharacteristics struct {
	MeanEstimation, VarianceEstimation float64
	Size                               int
}

func SamplePoolCharacteristics(sc []SampleCharacteristics) (pool SampleCharacteristics) {
	for j := range sc {
		sizeJ := float64(sc[j].Size)
		pool.MeanEstimation += sc[j].MeanEstimation * sizeJ
		pool.Size += sc[j].Size
		pool.VarianceEstimation += sc[j].VarianceEstimation * sizeJ * sizeJ
	}
	sizeTotal := float64(pool.Size)
	pool.MeanEstimation /= sizeTotal
	pool.VarianceEstimation /= sizeTotal * sizeTotal
	return
}

func SameSizePoolVariance(vars []float64) float64 {
	return SumFloat64Slice(vars) / (float64(len(vars)) * float64(len(vars)))
}

func LinearlyWeightedSameSizePoolMean(means []float64) (wm float64) {
	sumTerms := make([]float64, len(means))
	scale := float64(len(means) * (len(means) + 1) / 2)
	for i := range means {
		sumTerms[i] = float64(i+1) * means[i]
	}
	wm = SumFloat64Slice(sumTerms)
	wm /= scale
	return wm
}

func LinearlyWeightedSameSizePoolVariances(variances []float64) (wv float64) {
	sumTerms := make([]float64, len(variances))
	scale := float64(len(variances) * (len(variances) + 1) / 2)
	scale *= scale
	for i := range variances {
		sumTerms[i] = float64((i+1)*(i+1)) * variances[i]
	}
	wv = SumFloat64Slice(sumTerms)
	wv /= scale
	return wv
}
