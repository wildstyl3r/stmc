package extensions

import (
	"math"

	"github.com/wildstyl3r/lxgata"
	"github.com/wildstyl3r/stmc/internal/constants"
	"github.com/wildstyl3r/stmc/internal/model"
	"github.com/wildstyl3r/stmc/internal/utils"
)

const GlowDischargeDensityPPHCKey = "GlowDischargeDensityPPHC"

func GlowDischargeDensityPPHC(m *model.Model) ([]string, []any, error) {
	ionizations := utils.GriddedInterval{
		Values: m.GetMetrics(model.DataHubKeyType(NormalizedCollisionRateKey)).(map[lxgata.CollisionType][]float64)[lxgata.IONIZATION],
		Step:   m.XStep}
	// L = 2*H
	// n(x) = ( A / 2D ) * ([-x] / ([2d] + [2H]) ) * ( [2H]{d,x}p - [2d]{x,H}p - [2(d+H)]{d,H}m - [2(x+d)]{d,x}m + [2(x+H)]{x,H}m + [2x]{d,H}p )

	var fluxAtCathode float64
	if m.Parameters.CalculateCurrentDensity {
		fluxAtCathode = 1 //dummy flux just to find density shape
	} else {
		fluxAtCathode = m.Parameters.CathodeCurrentDensity / constants.ElementaryCharge
	}

	Lambda := m.Parameters.TubeRadius / 2.4
	powerFactor := 1. / Lambda
	ambipolarDiffusionCoefficient := m.Parameters.IonMobility * m.Parameters.SlowElectronTemperature
	scaleFactor := Lambda / (2 * ambipolarDiffusionCoefficient)

	expBrace := func(v float64) float64 { return math.Exp(v * powerFactor) }

	exp2H, exp2d := expBrace(2*m.Parameters.GapLength), expBrace(2*m.Parameters.CathodeFallLength) //GapLength in case of PPHC is already halved in NewModel()
	normDivisor := exp2d + exp2H

	minusTerms, plusTerms, density := make([]float64, len(ionizations.Values)), make([]float64, len(ionizations.Values)), make([]float64, len(ionizations.Values))

	dIndex := int(math.Floor(m.Parameters.CathodeFallLength / m.XStep))
	for i := dIndex; i < len(ionizations.Values); i++ {
		x := m.XStep * float64(i)
		minusTerms[i] = expBrace(-x) * ionizations.Interpolate(x)
		plusTerms[i] = expBrace(x) * ionizations.Interpolate(x)
	}

	thirdOuterSumTerm := exp2H * exp2d * utils.SumFloat64Slice(minusTerms[dIndex:])
	plusFullIntegral := utils.SumFloat64Slice(plusTerms[dIndex:])

	for i := dIndex; i < len(density); i++ {
		x := m.XStep * float64(i)
		exp2x := expBrace(2 * x)
		density[i] = max(0, fluxAtCathode*scaleFactor*expBrace(-x)/normDivisor*utils.SumFloat64Slice([]float64{
			exp2H * utils.SumFloat64Slice(plusTerms[dIndex:i]),
			-exp2d * utils.SumFloat64Slice(plusTerms[i:]),
			thirdOuterSumTerm,
			exp2x * exp2d * utils.SumFloat64Slice(minusTerms[dIndex:i]),
			exp2x * exp2H * utils.SumFloat64Slice(minusTerms[i:]),
			exp2x * plusFullIntegral,
		}))
	}
	return []string{GlowDischargeDensityPPHCKey}, []any{density}, nil
}
