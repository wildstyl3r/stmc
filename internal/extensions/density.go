package extensions

import (
	"math"
	"sort"

	"github.com/wildstyl3r/lxgata"
	"github.com/wildstyl3r/stmc/internal/model"
	"github.com/wildstyl3r/stmc/internal/utils"
)

const (
	GlowDischargeDensityKey           = "GlowDischargeDensity"
	SimplifiedGlowDischargeDensityKey = "SimplifiedGlowDischargeDensity"
)

func GlowDischargeDensity(m *model.Model) ([]string, []any, error) {
	ionizations := utils.GriddedInterval{
		Values: m.Get(model.DataHubKeyType(NormalizedCollisionRateKey)).(map[lxgata.CollisionType][]float64)[lxgata.IONIZATION],
		Step:   m.XStep}
	density := make([]float64, len(ionizations.Values))
	// n(x) = A * [-x] / (2D ([2L] - [2d])) * ([2L]B - [2x]E + [2d]F - [2(d+L)]G + [2(x+d)]K + [2(x+L)]M)
	// B = int+_d^x,   E = int+_d^L,   F = int+_x^L, G = int-_d^L,   K = int-_x^L,   M = int-_x^L

	Lambda := m.Parameters.CathodeRadius / 2.4
	powerFactor := 1. / Lambda

	expBrace := func(v float64) float64 { return math.Exp(v * powerFactor) }
	plusTypeIntegrable := func(v float64) float64 {
		return expBrace(v) * ionizations.Interpolate(v)
	}
	minusTypeIntegrable := func(v float64) float64 {
		return expBrace(-v) * ionizations.Interpolate(v)
	}

	evalPlusZetaAtDc, evalPlusZetaAtL := plusTypeIntegrable(m.Parameters.CathodeFallLength), plusTypeIntegrable(m.Parameters.GapLength)

	dcCeilIndex := int(math.Ceil(m.Parameters.CathodeFallLength / m.XStep))
	LFloorIndex := int(math.Floor(m.Parameters.GapLength / m.XStep))
	// summing from smaller to greater values
	accum := plusTypeIntegrable((float64(dcCeilIndex)+0.5)*m.XStep) * (float64(dcCeilIndex)*m.XStep - m.Parameters.CathodeFallLength)
	plusAtL := accum + plusTypeIntegrable(m.Parameters.GapLength)*(m.Parameters.GapLength-float64(LFloorIndex)*m.XStep)
	plusPrefixSums := make([]float64, len(ionizations.Values)+1)
	if evalPlusZetaAtDc < evalPlusZetaAtL {
		for i := dcCeilIndex + 1; i <= LFloorIndex; i++ {
			currentAddition := plusTypeIntegrable((float64(i)+0.5)*m.XStep) * m.XStep
			plusPrefixSums[i] = accum + currentAddition
			accum = plusPrefixSums[i]
			plusAtL += currentAddition
		}
	} else {
		for i := LFloorIndex; i > dcCeilIndex; i-- {
			currentAddition := plusTypeIntegrable((float64(i)+0.5)*m.XStep) * m.XStep
			plusPrefixSums[i] = accum + currentAddition
			accum = plusPrefixSums[i]
			plusAtL += currentAddition
		}
	}

	evalMinusZetaAtDc, evalMinusZetaAtL := minusTypeIntegrable(m.Parameters.CathodeFallLength), minusTypeIntegrable(m.Parameters.GapLength)
	accum = minusTypeIntegrable((float64(dcCeilIndex)+0.5)*m.XStep) * (float64(dcCeilIndex)*m.XStep - m.Parameters.CathodeFallLength)
	minusAtL := accum + minusTypeIntegrable(m.Parameters.GapLength)*(m.Parameters.GapLength-float64(LFloorIndex)*m.XStep)
	minusPrefixSums := make([]float64, len(ionizations.Values)+1)
	if evalMinusZetaAtDc < evalMinusZetaAtL {
		for i := dcCeilIndex + 1; i <= LFloorIndex; i++ {
			currentAddition := minusTypeIntegrable((float64(i)+0.5)*m.XStep) * m.XStep
			minusPrefixSums[i] = accum + currentAddition
			accum = minusPrefixSums[i]
			minusAtL += currentAddition
		}
	} else {
		for i := LFloorIndex; i > dcCeilIndex; i-- {
			currentAddition := minusTypeIntegrable((float64(i)+0.5)*m.XStep) * m.XStep
			minusPrefixSums[i] = accum + currentAddition
			accum = minusPrefixSums[i]
			minusAtL += currentAddition
		}
	}

	commonFactor1 := Lambda / (2 * m.Parameters.AmbipolarDiffusionCoefficient)
	commonDivisor := expBrace(2*m.Parameters.GapLength) - expBrace(2*m.Parameters.CathodeFallLength)
	for i := dcCeilIndex; i <= LFloorIndex; i++ {
		positiveSumTerms := []float64{
			expBrace(2*m.Parameters.GapLength) * plusPrefixSums[i],                                        //B
			expBrace(2*m.Parameters.CathodeFallLength) * (plusAtL - plusPrefixSums[i]),                    //F
			expBrace(2*(m.Parameters.CathodeFallLength+float64(i)*m.XStep)) * minusPrefixSums[i],          //K
			expBrace(2 * (float64(i)*m.XStep + m.Parameters.GapLength) * (minusAtL - minusPrefixSums[i])), //M
		}
		negativeSumTerms := []float64{
			expBrace(2*float64(i)*m.XStep) * plusAtL,                                       //E
			expBrace(2*(m.Parameters.CathodeFallLength+m.Parameters.GapLength)) * minusAtL, //G
		}
		sort.Float64s(positiveSumTerms)
		sort.Float64s(negativeSumTerms)

		commonFactor2 := expBrace(-float64(i) * m.XStep)
		sumOfIntegrals := 0.
		for j := range positiveSumTerms {
			sumOfIntegrals += positiveSumTerms[j]
		}
		for j := range negativeSumTerms {
			sumOfIntegrals -= negativeSumTerms[j]
		}
		density[i] = commonFactor1 * commonFactor2 * sumOfIntegrals / commonDivisor
	}
	return []string{GlowDischargeDensityKey}, []any{density}, nil
}

func SimplifiedGlowDischargeDensity(m *model.Model) ([]string, []any, error) {
	ionizations := m.Get(model.DataHubKeyType(NormalizedCollisionRateKey)).(map[lxgata.CollisionType][]float64)[lxgata.IONIZATION]
	density := make([]float64, len(ionizations))

	ionizationsGI := utils.GriddedInterval{Values: append(ionizations, ionizations[len(ionizations)-1], ionizations[len(ionizations)-1]), Step: m.XStep}

	preC1, err := utils.VariableLimitDoubleIntegration(
		m.Parameters.CathodeFallLength, m.Parameters.GapLength, m.XStep,
		0, func(x float64) float64 { return x },
		ionizationsGI)
	if err != nil {
		return nil, nil, err
	}
	C1 := -preC1 / (m.Parameters.GapLength - m.Parameters.CathodeFallLength)

	preC2, err := utils.VariableLimitDoubleIntegration(
		0, m.Parameters.GapLength, m.XStep,
		0, func(x float64) float64 { return x },
		ionizationsGI)
	if err != nil {
		return nil, nil, err
	}
	C2 := -(preC2 + C1*m.Parameters.GapLength)
	{
		sum := 0.
		accum := 0.
		for j := range ionizations {
			sum += accum
			accum += ionizations[j] * m.XStep
			density[j] = -(sum*m.XStep + C1*float64(j)*m.XStep + C2)
		}
	}
	return []string{SimplifiedGlowDischargeDensityKey}, []any{density}, nil
}
