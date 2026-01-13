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

type WeightedAggregation struct {
	sum              KahanSummable
	squaresSum       KahanSummable
	weightSum        KahanSummable
	weightSquaresSum KahanSummable
}

func (a *WeightedAggregation) Update(values, weights []float64) {
	if weights == nil {
		for i := range values {
			a.sum.Add(values[i])
			a.weightSum.Add(1)
			a.squaresSum.Add(values[i] * values[i])
			a.weightSquaresSum.Add(1)
		}
	} else {
		for i := range values {
			a.sum.Add(values[i] * weights[i])
			a.weightSum.Add(weights[i])
			a.squaresSum.Add(values[i] * values[i] * weights[i])
			a.weightSquaresSum.Add(weights[i] * weights[i])
		}
	}
	// a.sum.Add(SumFloat64Slice(values, true))
	// a.WeightSum.Add(SumFloat64Slice(weights, true))
}

func (a *WeightedAggregation) UpdateK(values, weights []KahanSummable) {
	if weights == nil {
		for i := range values {
			v := values[i].Val()
			a.sum.Add(v)
			a.weightSum.Add(1)
			a.squaresSum.Add(v * v)
			a.weightSquaresSum.Add(1)
		}
	} else {
		for i := range values {
			v, w := values[i].Val(), weights[i].Val()
			a.sum.Add(v * w)
			a.weightSum.Add(w)
			a.squaresSum.Add(v * v * w)
			a.weightSquaresSum.Add(w * w)
		}
	}
	// a.sum.Add(SumFloat64Slice(values, true))
	// a.WeightSum.Add(SumFloat64Slice(weights, true))
}
func (a *WeightedAggregation) UpdateIK(values []int, weights []KahanSummable) {
	if weights == nil {
		for i := range values {
			v := float64(values[i])
			a.sum.Add(v)
			if math.IsNaN(a.sum.val) {
				println("a.sum  is nan")
			}
			a.weightSum.Add(1)
			a.squaresSum.Add(v * v)
			a.weightSquaresSum.Add(1)
		}
	} else {
		for i := range values {
			v, w := float64(values[i]), weights[i].Val()
			a.sum.Add(v * w)
			if math.IsNaN(a.sum.val) {
				println("a.sum  is nan")
			}
			a.weightSum.Add(w)
			a.squaresSum.Add(v * v * w)
			a.weightSquaresSum.Add(w * w)
		}
	}

	// a.sum.Add(SumFloat64Slice(values, true))
	// a.WeightSum.Add(SumFloat64Slice(weights, true))
}

func (a *WeightedAggregation) IncrementBy1() {
	a.squaresSum.Add(2*a.sum.Val() + a.weightSum.Val())
	a.sum.Add(a.weightSum.Val())
}

func (a *WeightedAggregation) Mean() float64 {
	return a.sum.Val() / a.weightSum.Val()
}

func (a *WeightedAggregation) meanAndVariance() (float64, float64) {
	if a.sum.val == 0 {
		return 0, 0
	}
	return a.Mean(), (a.weightSum.Val()*a.squaresSum.Val() - a.sum.Val()*a.sum.Val()) / (a.weightSum.Val()*a.weightSum.Val() - a.weightSquaresSum.Val())
}

func (a *WeightedAggregation) MeanWithErrorMargin(confidenceP float64) (float64, float64) {
	m, v := a.meanAndVariance()
	if v == 0 || math.IsNaN(v) || math.IsInf(v, 0) {
		return m, 0
	}
	return m, distuv.UnitNormal.Quantile(1.-(1.-confidenceP)*0.5) * math.Sqrt(v/a.EffectiveSampleSize()) //EstimateMargin(confidenceP, v, a.EffectiveSampleSize())
}

func (a *WeightedAggregation) EffectiveSampleSize() float64 {
	ws := a.weightSum.Val()
	return ws * ws / a.weightSquaresSum.Val()
}

type WeightedDistribution struct {
	Cells []KahanSummable
	Mass  KahanSummable
}

func NewWeightedDistribution(cellNumber int) WeightedDistribution {
	return WeightedDistribution{
		Cells: make([]KahanSummable, cellNumber),
	}
}

type WeightedDistributionUpdate struct {
	cell  int
	value float64
}

func (a *WeightedDistribution) Update(update []WeightedDistributionUpdate, weights []float64) {
	for element := range update {
		a.Cells[update[element].cell].Add(update[element].value * weights[element])
	}
	a.Mass.Add(SumFloat64Slice(weights, true))
}

func (a *WeightedDistribution) Mean(cell int) float64 {
	if a.Mass.Val() == 0 {
		return 0
	}
	return a.Cells[cell].Val() / a.Mass.Val()
}

type AggregatedDistribution struct {
	Cells []AggregationF
}

func NewAggregatedDistribution(cellNumber int) AggregatedDistribution {
	return AggregatedDistribution{
		Cells: make([]AggregationF, cellNumber),
	}
}

func (d *AggregatedDistribution) Update(values []WeightedDistribution) {
	for cell := range d.Cells {
		cellValues := make([]float64, len(values))
		for sample := range values {
			cellValues[sample] = values[sample].Mean(cell)
		}
		d.Cells[cell].Update(cellValues)
	}
}

type AggregationF struct {
	sum       KahanSummable
	squareSum KahanSummable
	Number    int

	// devSquareSum KahanSummable
}

func (a *AggregationF) Append(piece KahanSummable) {
	v := piece.Val()
	a.sum.Add(v)
	a.Number++
	a.squareSum.Add(v * v)
	// a.devSquareSum.Add()
}
func (a *AggregationF) AppendF(piece float64) {
	a.sum.Add(piece)
	a.Number++
	a.squareSum.Add(piece * piece)
	// a.devSquareSum.Add()
}

func (a *AggregationF) MeanWithErrorMargin(confidenceP float64) (float64, float64) {
	m, v := a.meanAndVariance()
	return m, EstimateMargin(confidenceP, v, float64(a.Number))
}

func (a *AggregationF) Update(values []float64) {
	a.sum.Add(SumFloat64Slice(values, true))
	a.Number += len(values)
	// mean := a.sum.Val() / float64(a.Number)
	// dsSumTerms := make([]float64, len(values))
	for i := range values {
		// dev := mean - values[i]
		// dsSumTerms[i] = dev * dev
		a.squareSum.Add(values[i] * values[i])
	}

	// a.devSquareSum.Add(SumFloat64Slice(dsSumTerms, true))

}

func (a *AggregationF) UpdateK(values []KahanSummable) {
	terms := make([]float64, len(values))
	for i := range terms {
		terms[i] = values[i].Val()
	}
	a.sum.Add(SumFloat64Slice(terms, true))
	a.Number += len(values)
	// mean := a.sum.Val() / float64(a.Number)
	// dsSumTerms := make([]float64, len(values))
	for i := range values {
		// dev := mean - terms[i]
		// dsSumTerms[i] = dev * dev
		v := values[i].Val()
		a.squareSum.Add(v * v)
	}

	// a.devSquareSum.Add(SumFloat64Slice(dsSumTerms, true))
}

func (a *AggregationF) Mean() float64 {
	return a.sum.Val() / float64(a.Number)
}

func (a *AggregationF) meanAndVariance() (float64, float64) {
	//variance estimation may be slightly worse than in non-aggregated case, but we should be OK with that
	m := a.Mean()
	return m, a.squareSum.Val()/float64(a.Number) - m*m //a.Mean(), a.devSquareSum.Val() / float64(a.Number-1)
}

type Aggregation struct {
	sum, Number int
	squareSum   int

	// DevSquareSum KahanSummable
}

func (a *Aggregation) MeanWithErrorMargin(confidenceP float64) (float64, float64) {
	m, v := a.meanAndVariance()
	if v == 0 || math.IsNaN(v) || math.IsInf(v, 0) {
		return m, 0
	}
	return m, distuv.UnitNormal.Quantile(1.-(1.-confidenceP)*0.5) * math.Sqrt(v/float64(a.Number)) //EstimateMargin(confidenceP, v, a.EffectiveSampleSize())
}

func (a *Aggregation) IncrementBy1() {
	a.sum += a.Number
}

func (a *Aggregation) Update(values []int) {
	a.sum += SumIntSlice(values)
	a.Number += len(values)
	// mean := float64(a.sum) / float64(a.Number)
	// dsSumTerms := make([]float64, len(values))
	for i := range values {
		a.squareSum += values[i] * values[i]
	}
	// a.DevSquareSum.Add(SumFloat64Slice(dsSumTerms, true))
}

func (a *Aggregation) Mean() float64 {
	return float64(a.sum) / float64(a.Number)
}

func (a *Aggregation) meanAndVariance() (float64, float64) {
	m := a.Mean()
	return m, float64(a.squareSum)/float64(a.Number) - m*m
}

func Mean[T Number](s []T) (mean float64) {
	floated := make([]float64, len(s))
	for i := range s {
		floated[i] = float64(s[i])
	}
	mean = SumFloat64Slice(floated, true) / float64(len(s))
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
	return n * SumFloat64Slice(sumTerms, true) / ((n - 1) * (n - 2))
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
		raw[r-1] = SumFloat64Slice(sumTerms, true) / n
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
	return SumFloat64Slice(sumTerms, true)
}

func StudentedMarginFromData(confidenceP float64, data []float64) float64 {
	dist := distuv.StudentsT{Nu: float64(len(data) - 1), Sigma: 1}
	quantileFactor := dist.Quantile(1 - (1-confidenceP)/2) // make it two-sided
	return quantileFactor * math.Sqrt(Variance(data, true)/float64(len(data)))
}

func EstimateMargin(confidenceP, variance, sampleSize float64) float64 {
	var quantile float64
	twoSidedP := 1. - (1.-confidenceP)*0.5
	if sampleSize > 33 {
		quantile = distuv.UnitNormal.Quantile(twoSidedP)
	} else if sampleSize > 1 {
		quantile = distuv.StudentsT{Nu: sampleSize - 1, Sigma: 1}.Quantile(twoSidedP)
	} else {
		return 5 * math.Sqrt(variance)
	}
	se := math.Sqrt(variance / sampleSize)
	return quantile * se
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
	return SumFloat64Slice(vars, true) / (float64(len(vars)) * float64(len(vars)))
}

func LinearlyWeightedSameSizePoolMean(means []float64) (wm float64) {
	sumTerms := make([]float64, len(means))
	scale := float64(len(means) * (len(means) + 1) / 2)
	for i := range means {
		sumTerms[i] = float64(i+1) * means[i]
	}
	wm = SumFloat64Slice(sumTerms, true)
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
	wv = SumFloat64Slice(sumTerms, true)
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
	b1 = SumFloat64Slice(numeratorTerms, true) / SumFloat64Slice(denominatorTerms, true)
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
	denominator := SumFloat64Slice(denominatorTerms, true)
	b1 = SumFloat64Slice(numeratorTerms, true) / denominator
	b1Variance = SumFloat64Slice(b1VarianceTerms, true) / (denominator * denominator)
	b0 = meanY - b1*meanX
	b0Variance = b1Variance*(meanX*meanX) + SumFloat64Slice(yVariance, true)/float64(n*n)
	return
}

// func QuadraticRegressionMSE(x, y []float64) (b2, b1, b0 float64) { //y ~ b_2*x + b_1*x + b_0 + N(0,?)

// }
