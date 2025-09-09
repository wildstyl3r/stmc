package extensions

import (
	"math"
	"sort"

	"github.com/wildstyl3r/lxgata"
	"github.com/wildstyl3r/stmc/internal/datahub"
	"github.com/wildstyl3r/stmc/internal/model"
	"github.com/wildstyl3r/stmc/internal/utils"
)

const (
	GlowDischargeDensityKey           = "GlowDischargeDensity"
	SimplifiedGlowDischargeDensityKey = "SimplifiedGlowDischargeDensity"
)

func GlowDischargeDensity(model *model.Model) ([]string, []any, error) {
	ionizations := utils.GriddedInterval{
		Values: datahub.Get(datahub.KeyType(lxgata.IONIZATION), model).([]float64),
		Step:   model.XStep}
	density := make([]float64, len(ionizations.Values))
	// n(x) = A * [-x] / (2D ([2L] - [2d])) * ([2L]B - [2x]E + [2d]F - [2(d+L)]G + [2(x+d)]K + [2(x+L)]M)
	// B = int+_d^x,   E = int+_d^L,   F = int+_x^L, G = int-_d^L,   K = int-_x^L,   M = int-_x^L

	Lambda := model.Parameters.CathodeRadius / 2.4
	powerFactor := 1. / Lambda

	expBrace := func(v float64) float64 { return math.Exp(v * powerFactor) }
	plusTypeIntegrable := func(v float64) float64 {
		return expBrace(v) * ionizations.Interpolate(v)
	}
	minusTypeIntegrable := func(v float64) float64 {
		return expBrace(-v) * ionizations.Interpolate(v)
	}

	evalPlusZetaAtDc, evalPlusZetaAtL := plusTypeIntegrable(model.Parameters.CathodeFallLength), plusTypeIntegrable(model.Parameters.GapLength)

	dcCeilIndex := int(math.Ceil(model.Parameters.CathodeFallLength / model.XStep))
	LFloorIndex := int(math.Floor(model.Parameters.GapLength / model.XStep))
	// summing from smaller to greater values
	accum := plusTypeIntegrable((float64(dcCeilIndex)+0.5)*model.XStep) * (float64(dcCeilIndex)*model.XStep - model.Parameters.CathodeFallLength)
	plusAtL := accum + plusTypeIntegrable(model.Parameters.GapLength)*(model.Parameters.GapLength-float64(LFloorIndex)*model.XStep)
	plusPrefixSums := make([]float64, len(ionizations.Values)+1)
	if evalPlusZetaAtDc < evalPlusZetaAtL {
		for i := dcCeilIndex + 1; i <= LFloorIndex; i++ {
			currentAddition := plusTypeIntegrable((float64(i)+0.5)*model.XStep) * model.XStep
			plusPrefixSums[i] = accum + currentAddition
			accum = plusPrefixSums[i]
			plusAtL += currentAddition
		}
	} else {
		for i := LFloorIndex; i > dcCeilIndex; i-- {
			currentAddition := plusTypeIntegrable((float64(i)+0.5)*model.XStep) * model.XStep
			plusPrefixSums[i] = accum + currentAddition
			accum = plusPrefixSums[i]
			plusAtL += currentAddition
		}
	}

	evalMinusZetaAtDc, evalMinusZetaAtL := minusTypeIntegrable(model.Parameters.CathodeFallLength), minusTypeIntegrable(model.Parameters.GapLength)
	accum = minusTypeIntegrable((float64(dcCeilIndex)+0.5)*model.XStep) * (float64(dcCeilIndex)*model.XStep - model.Parameters.CathodeFallLength)
	minusAtL := accum + minusTypeIntegrable(model.Parameters.GapLength)*(model.Parameters.GapLength-float64(LFloorIndex)*model.XStep)
	minusPrefixSums := make([]float64, len(ionizations.Values)+1)
	if evalMinusZetaAtDc < evalMinusZetaAtL {
		for i := dcCeilIndex + 1; i <= LFloorIndex; i++ {
			currentAddition := minusTypeIntegrable((float64(i)+0.5)*model.XStep) * model.XStep
			minusPrefixSums[i] = accum + currentAddition
			accum = minusPrefixSums[i]
			minusAtL += currentAddition
		}
	} else {
		for i := LFloorIndex; i > dcCeilIndex; i-- {
			currentAddition := minusTypeIntegrable((float64(i)+0.5)*model.XStep) * model.XStep
			minusPrefixSums[i] = accum + currentAddition
			accum = minusPrefixSums[i]
			minusAtL += currentAddition
		}
	}

	commonFactor1 := Lambda / (2 * model.Parameters.AmbipolarDiffusionCoefficient)
	commonDivisor := expBrace(2*model.Parameters.GapLength) - expBrace(2*model.Parameters.CathodeFallLength)
	for i := dcCeilIndex; i <= LFloorIndex; i++ {
		positiveSumTerms := []float64{
			expBrace(2*model.Parameters.GapLength) * plusPrefixSums[i],                                            //B
			expBrace(2*model.Parameters.CathodeFallLength) * (plusAtL - plusPrefixSums[i]),                        //F
			expBrace(2*(model.Parameters.CathodeFallLength+float64(i)*model.XStep)) * minusPrefixSums[i],          //K
			expBrace(2 * (float64(i)*model.XStep + model.Parameters.GapLength) * (minusAtL - minusPrefixSums[i])), //M
		}
		negativeSumTerms := []float64{
			expBrace(2*float64(i)*model.XStep) * plusAtL,                                           //E
			expBrace(2*(model.Parameters.CathodeFallLength+model.Parameters.GapLength)) * minusAtL, //G
		}
		sort.Float64s(positiveSumTerms)
		sort.Float64s(negativeSumTerms)

		commonFactor2 := expBrace(-float64(i) * model.XStep)
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

func SimplifiedGlowDischargeDensity(model *model.Model) ([]string, []any, error) {
	ionizations := datahub.Get(datahub.KeyType(lxgata.IONIZATION), model).([]float64)
	density := make([]float64, len(ionizations))

	ionizationsGI := utils.GriddedInterval{Values: append(ionizations, ionizations[len(ionizations)-1]), Step: model.XStep}

	preC1, err := utils.VariableLimitDoubleIntegration(
		model.Parameters.CathodeFallLength, model.Parameters.GapLength, model.XStep,
		0, func(x float64) float64 { return x },
		ionizationsGI)
	if err != nil {
		return nil, nil, err
	}
	C1 := -preC1 / (model.Parameters.GapLength - model.Parameters.CathodeFallLength)

	preC2, err := utils.VariableLimitDoubleIntegration(
		0, model.Parameters.GapLength, model.XStep,
		0, func(x float64) float64 { return x },
		ionizationsGI)
	if err != nil {
		return nil, nil, err
	}
	C2 := -(preC2 + C1*model.Parameters.GapLength)
	{
		sum := 0.
		accum := 0.
		for j := range ionizations {
			sum += accum
			accum += ionizations[j] * model.XStep
			density[j] = -(sum*model.XStep + C1*float64(j)*model.XStep + C2)
		}
	}
	return []string{SimplifiedGlowDischargeDensityKey}, []any{density}, nil
}
