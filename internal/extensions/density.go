package extensions

import (
	"math"

	"github.com/wildstyl3r/lxgata"
	"github.com/wildstyl3r/stmc/internal/datahub"
	"github.com/wildstyl3r/stmc/internal/model"
	"github.com/wildstyl3r/stmc/internal/utils"
)

const GlowDischargeDensityKey = "GlowDischargeDensity"
const SimplifiedGlowDischargeDensityKey = "SimplifiedGlowDischargeDensity"

func GlowDischargeDensity(model *model.Model) ([]string, []any) {
	ionizations := datahub.Get(NormalizedCollisionRateKey, model).(map[lxgata.CollisionType][]float64)[lxgata.IONIZATION]
	density := make([]float64, len(ionizations))
	Lambda := model.Parameters.CathodeRadius / 2.4
	powerFactor := 1. / Lambda
	integralFactor := 1. //2. * Lambda / model.Parameters.AmbipolarDiffusionCoefficient

	a := math.Exp(2 * model.Parameters.GapLength * powerFactor)
	e := math.Exp(2 * model.Parameters.CathodeFallLength * powerFactor)
	g := math.Exp(2 * (model.Parameters.CathodeFallLength + model.Parameters.GapLength) * powerFactor)
	scaler := 1. / (e - a)

	prefSumMinus := make([]float64, len(ionizations))
	prefSumPlus := make([]float64, len(ionizations))
	for i := range len(ionizations) {
		x := float64(i) * model.XStep
		prefSumMinus[i] = math.Exp(-x*powerFactor) * ionizations[i]
		prefSumPlus[i] = math.Exp(x*powerFactor) * ionizations[i]
		if i > 0 {
			prefSumMinus[i] += prefSumMinus[i-1]
			prefSumPlus[i] += prefSumPlus[i-1]
		}
	}
	for i := range len(ionizations) {
		prefSumMinus[i] *= model.XStep
		prefSumPlus[i] *= model.XStep
	}
	L := len(prefSumPlus) - 1
	d := min(int(model.Parameters.CathodeFallLength/model.XStep), L-1)
	IPlus := prefSumPlus[L] - prefSumPlus[d]
	IMinus := prefSumMinus[L] - prefSumMinus[d]
	for i := 0; i < len(ionizations); i++ {
		x := float64(i) * model.XStep
		density[i] = scaler * math.Exp(-x*powerFactor) * integralFactor * (-a*(prefSumPlus[i]-prefSumPlus[d]) +
			math.Exp(2*x*powerFactor)*IPlus -
			e*(prefSumPlus[L]-prefSumPlus[i]) +
			g*IMinus -
			math.Exp(2*(model.Parameters.GapLength+x)*powerFactor)*(prefSumMinus[L]-prefSumMinus[i]) -
			math.Exp(2*(model.Parameters.CathodeFallLength+x)*powerFactor)*(prefSumMinus[i]-prefSumMinus[d]))
	}
	return []string{GlowDischargeDensityKey}, []any{density}
}

func SimplifiedGlowDischargeDensity(model *model.Model) ([]string, []any) {
	ionizations := datahub.Get(datahub.KeyType(lxgata.IONIZATION), model).([]float64)
	density := make([]float64, len(ionizations))

	dcIndex := min(int(model.Parameters.CathodeFallLength/model.XStep), model.NumCells-2)
	exactSnAtDc := ionizations[dcIndex] + (model.Parameters.CathodeFallLength-float64(dcIndex)*model.XStep)/
		model.XStep*(ionizations[dcIndex+1]-ionizations[dcIndex])
	C1 := 0.5 * (exactSnAtDc + ionizations[dcIndex]) * model.XStep
	{
		accum := utils.TableIntegrate(ionizations[:dcIndex], nil, model.XStep)
		for i := dcIndex; i < len(ionizations); i++ {
			C1 += accum
			accum += ionizations[i] * model.XStep
		}
		C1 *= -model.XStep
		C1 /= (model.Parameters.GapLength - model.Parameters.CathodeFallLength)
	}
	C2 := 0.
	{
		accum := 0.
		for i := range ionizations {
			C2 -= accum
			accum += ionizations[i] * model.XStep //TableIntegrate(datahub.GetKey(datahub.KeyType(lxgata.IONIZATION))[:i], nil, model.XStep)
		}
		C2 *= model.XStep
		C2 -= C1 * model.Parameters.GapLength
	}
	{
		sum := 0.
		accum := 0.
		for j := range ionizations {
			sum += accum
			accum += ionizations[j] * model.XStep
			if j > dcIndex || true {
				density[j] = -(sum*model.XStep + C1*float64(j)*model.XStep + C2)
			}
		}
	}
	return []string{SimplifiedGlowDischargeDensityKey}, []any{density}
}
