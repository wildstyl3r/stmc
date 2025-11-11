package utils

import (
	"math"
	"slices"

	"golang.org/x/exp/constraints"
	"gonum.org/v1/gonum/stat/distuv"
)

type Number interface {
	constraints.Float | constraints.Integer
}

type Aggregation struct {
	Sum, Number int

	DevSquareSum float64
}

func (a *Aggregation) Update(values []int) {
	a.Sum += SumIntSlice(values)
	a.Number += len(values)
	mean := float64(a.Sum) / float64(a.Number)
	dsSumTerms := make([]float64, len(values))
	for i := range values {
		dev := mean - float64(values[i])
		dsSumTerms[i] = dev * dev
	}
	a.DevSquareSum += SumFloat64Slice(dsSumTerms)
}

func (a *Aggregation) Mean() float64 {
	return float64(a.Sum) / float64(a.Number)
}

func (a *Aggregation) MeanAndVariance() (float64, float64) {
	//variance estimation may be slightly worse than in non-aggregated case, but we should be OK with that
	return a.Mean(), a.DevSquareSum / float64(a.Number-1)
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

func NormalMargin(confidenceP, variance float64, sampleSize int) float64 {
	return distuv.UnitNormal.Quantile(1.-(1.-confidenceP)*0.5) * math.Sqrt(variance/float64(sampleSize))
}

func WilsonHilfertyCoefficient[T constraints.Integer](confidenceP float64, k T, lower bool) float64 {
	kf := float64(k)
	if lower {
		return float64(kf) * math.Pow(1.-1./(9*kf)-
			distuv.UnitNormal.Quantile(1.-(1.-confidenceP)*0.5)/(3*math.Sqrt(kf)), 3)
	} else {
		return float64(kf) * math.Pow(1.-1./(9*kf)+
			distuv.UnitNormal.Quantile(1.-(1.-confidenceP)*0.5)/(3*math.Sqrt(kf)), 3)
	}
}

func PoissonMargin[T constraints.Integer](confidenceP float64, observations []T) (left, right float64) {
	alpha := 1. - confidenceP
	k := SumIntSlice(observations)
	kf := float64(k)
	n := float64(len(observations))
	if k == 0 {
		return 0, -math.Log(alpha) / n //WilsonHilfertyCoefficient(twoSidedConfP, 1, false) / n
	}
	return distuv.Gamma{Alpha: kf, Beta: 1}.Quantile(0.5*alpha) / n, distuv.Gamma{Alpha: kf + 1, Beta: 1}.Quantile(1-0.5*alpha) / n //WilsonHilfertyCoefficient(twoSidedConfP, k, true) / n, WilsonHilfertyCoefficient(twoSidedConfP, k+1, false) / n
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

func RelativeExcessLeap(sample []float64) float64 {
	diff := 0.
	for i := range sample {
		if i > 0 {
			diff += math.Abs(sample[i] - sample[i-1])
		}
	}
	return diff / slices.Max(sample)
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
func EstimateMargin(confidenceP, variance float64, sampleSize int) float64 {
	if sampleSize > 33 {
		return NormalMargin(confidenceP, variance, sampleSize)
	} else {
		return StudentedMargin(confidenceP, variance, sampleSize)
	}
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

func LinearRegressionMSE(x, y []float64) (b1, b0 float64) { // y ~ b_1*x + b_0 + N(0,?)
	meanX, meanY := Mean(x), Mean(y)
	n := min(len(x), len(y))
	numeratorTerms, denominatorTerms := make([]float64, n), make([]float64, n)
	for i := range n {
		numeratorTerms[i] = (x[i] - meanX) * (y[i] - meanY)
		denominatorTerms[i] = (x[i] - meanX) * (x[i] - meanX)
	}
	b1 = SumFloat64Slice(numeratorTerms) / SumFloat64Slice(denominatorTerms)
	b0 = meanY - b1*meanX
	return
}

func LinearRegressionMSEWithVariance(x, y, yVariance []float64) (b1, b0, b1Variance, b0Variance float64) { // y ~ b_1*x + b_0 + N(0,?)
	meanX, meanY := Mean(x), Mean(y)
	n := min(len(x), len(y))
	numeratorTerms, denominatorTerms, b1VarianceTerms := make([]float64, n), make([]float64, n), make([]float64, n)
	for i := range n {
		deltaX := x[i] - meanX
		numeratorTerms[i] = deltaX * (y[i] - meanY)
		denominatorTerms[i] = deltaX * deltaX
		b1VarianceTerms[i] = yVariance[i] * (deltaX * deltaX)
	}
	denominator := SumFloat64Slice(denominatorTerms)
	b1 = SumFloat64Slice(numeratorTerms) / denominator
	b1Variance = SumFloat64Slice(b1VarianceTerms) / (denominator * denominator)
	b0 = meanY - b1*meanX
	b0Variance = b1Variance*(meanX*meanX) + SumFloat64Slice(yVariance)/float64(n*n)
	return
}

// func QuadraticRegressionMSE(x, y []float64) (b2, b1, b0 float64) { //y ~ b_2*x + b_1*x + b_0 + N(0,?)

// }
