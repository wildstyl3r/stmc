package extensions

import (
	"math"

	"github.com/wildstyl3r/lxgata"
	"github.com/wildstyl3r/stmc/internal/model"
	"github.com/wildstyl3r/stmc/internal/utils"
)

func DensityPerCathodeElectronFluxPPHC(m *model.Model) ([]string, []any, error) {
	ionizations := m.GetMetrics(model.DataHubKeyType(SingleElectronCollisionRateKey)).(map[lxgata.CollisionType]utils.GriddedInterval)[lxgata.IONIZATION]

	density := utils.GriddedInterval{
		Step:   m.XStep,
		Offset: ionizations.Offset,
		Values: make([]float64, m.NumCells),
	}
	ambipolarDiffusionCoefficient := m.GetMetrics(AmbipolarDiffusionCoefficientKey).(float64)

	// L = 2*H
	H := m.Parameters.SimulationLength()
	if m.Parameters.NoDiffusionLoss {
		scaledIonizations := ionizations.Scale(1. / ambipolarDiffusionCoefficient)
		C1, err := scaledIonizations.TrapezoidIntegration(m.Parameters.CathodeFallLength, H)
		if err != nil {
			return nil, nil, err
		}
		C2 := -C1 * m.Parameters.CathodeFallLength
		for i := range density.Values {
			x := float64(i) * m.XStep
			integral, err := scaledIonizations.VariableLimitDoubleIntegration(m.Parameters.CathodeFallLength, x, m.XStep, m.Parameters.CathodeFallLength, func(f float64) float64 { return f })
			if err != nil {
				return nil, nil, err
			}
			density.Values[i] = max(0, -integral+C1*x+C2)
		}
	} else {
		// L = 2*H
		// n(x) = ( A / 2D ) * ([-x] / ([2d] + [2H]) ) * ( [2H]{d,x}p - [2d]{x,H}p - [2(d+H)]{d,H}m - [2(x+d)]{d,x}m + [2(x+H)]{x,H}m + [2x]{d,H}p )

		Lambda := m.GetMetrics(CharacteristicDiffusionScaleKey).(float64)
		powerFactor := 1. / Lambda
		scaleFactor := Lambda / (2 * ambipolarDiffusionCoefficient)

		expBrace := func(v float64) float64 {
			return math.Exp(v * powerFactor)
		}

		exp2H, exp2d := expBrace(2*H), expBrace(2*m.Parameters.CathodeFallLength)
		normDivisor := exp2d + exp2H

		negExpBrace := func(x float64) float64 { return expBrace(-x) }
		minusTerms := ionizations.MulPointwise(negExpBrace)
		plusTerms := ionizations.MulPointwise(expBrace)

		fullMIntegral, err := minusTerms.TrapezoidIntegration(m.Parameters.CathodeFallLength, H)
		if err != nil {
			return nil, nil, err
		}

		fullPIntegral, err := plusTerms.TrapezoidIntegration(m.Parameters.CathodeFallLength, H)
		if err != nil {
			return nil, nil, err
		}

		thirdOuterSumTerm := exp2H * exp2d * fullMIntegral

		for i := range density.Values {
			x := m.XStep * float64(i)
			exp2x := expBrace(2 * x)
			firstTerm, err := plusTerms.TrapezoidIntegration(m.Parameters.CathodeFallLength, x)
			if err != nil {
				return nil, nil, err
			}
			secondTerm, err := plusTerms.TrapezoidIntegration(x, H)
			if err != nil {
				return nil, nil, err
			}
			fourthTerm, err := minusTerms.TrapezoidIntegration(m.Parameters.CathodeFallLength, x)
			if err != nil {
				return nil, nil, err
			}
			fifthTerm, err := minusTerms.TrapezoidIntegration(x, H)
			if err != nil {
				return nil, nil, err
			}
			density.Values[i] = max(0, scaleFactor*expBrace(-x)/normDivisor*utils.SumFloat64Slice([]float64{
				exp2H * firstTerm,
				-exp2d * secondTerm,
				-thirdOuterSumTerm,
				-exp2x * exp2d * fourthTerm,
				exp2x * exp2H * fifthTerm,
				exp2x * fullPIntegral,
			}))
		}
	}

	return []string{DensityPerCathodeElectronFluxKey}, []any{density}, nil
}
