package extensions

import (
	"math"

	"github.com/wildstyl3r/lxgata"
	"github.com/wildstyl3r/stmc/internal/constants"
	"github.com/wildstyl3r/stmc/internal/model"
	"github.com/wildstyl3r/stmc/internal/utils"
)

const (
	GlowDischargeDensityKey                         = "GlowDischargeDensity"
	CathodeElectronFluxKey                          = "CathodeElectronFlux"
	CathodeElectronFluxMarginKey                    = "CathodeElectronFluxMargin"
	CathodeIonFluxKey                               = "CathodeIonFlux"
	CathodeIonFluxMarginKey                         = "CathodeIonFluxMargin"
	FirstDensityPeakIndexKey                        = "FirstDensityPeakIndex"
	SourceIntegralPerCathodeElectronFluxKey         = "SingleElectronSourceIntegral"
	SourceIntegralPerCathodeElectronFluxVarianceKey = "SingleElectronSourceIntegralMargin"
	DensityPerCathodeElectronFluxKey                = "DensityPerCathodeElectronFluxKey"
)

func DensityPerCathodeElectronFlux(m *model.Model) ([]string, []any, error) {
	if m.Parameters.ParallelPlaneHollowCathode {
		return DensityPerCathodeElectronFluxPPHC(m)
	}
	ionizations := m.GetMetrics(model.DataHubKeyType(SingleElectronCollisionRateKey)).(map[lxgata.CollisionType]utils.GriddedInterval)[lxgata.IONIZATION]

	ambipolarDiffusionCoefficient := m.GetMetrics(AmbipolarDiffusionCoefficientKey).(float64)

	density := utils.GriddedInterval{
		Step:   m.XStep,
		Offset: ionizations.Offset,
		Values: make([]float64, m.NumCells),
	}

	if m.Parameters.DoNotAccountDiffusionLoss {
		// n(x) = -int{int{S(h),d,chi},d, x} + c1*x + c2
		scaledIonizations := ionizations.Scale(1. / ambipolarDiffusionCoefficient)
		C1, err := scaledIonizations.VariableLimitDoubleIntegration(m.Parameters.CathodeFallLength, m.Parameters.GapLength, m.XStep, 0, func(x float64) float64 { return x })
		if err != nil {
			return nil, nil, err
		}
		C1 /= (m.Parameters.GapLength - m.Parameters.CathodeFallLength)
		C2 := -C1 * m.Parameters.CathodeFallLength
		for i := range density.Values {
			x := float64(i) * m.XStep
			integral, err := scaledIonizations.VariableLimitDoubleIntegration(m.Parameters.CathodeFallLength, x, m.XStep, m.Parameters.CathodeFallLength, func(f float64) float64 { return f })
			if err != nil {
				return nil, nil, err
			}
			density.Values[i] = -integral + C1*x + C2
		}
	} else {
		// n(x) = ( A / 2D ) * ([-x] / ([2d] - [2L]) ) * ( [2(L+d)]{d,L}m + [2x]{d,L}p - [2L]{d,x}p - [2d]{x,L}p - [2(x+d)]{d,x}m - [2(x+L)]{x,L}m )
		// n(x) = scale * [-x] / norm  * ( [2(L+d)]fullMinus + [2x]fullPlus - [2L]{d,x}p - [2d]{x,L}p - [2(x+d)]{d,x}m - [2(x+L)]{x,L}m )
		Lambda := m.GetMetrics(CharacteristicDiffusionScaleKey).(float64)

		powerFactor := 1. / Lambda
		scaleFactor := Lambda / (2 * ambipolarDiffusionCoefficient)

		expBrace := func(v float64) float64 { return math.Exp(v * powerFactor) }
		exp2L, exp2d := expBrace(2*m.Parameters.GapLength), expBrace(2*m.Parameters.CathodeFallLength)
		normDivisor := exp2d - exp2L

		negExpBrace := func(x float64) float64 { return expBrace(-x) }
		minusTerms := ionizations.MulPointwise(negExpBrace)
		plusTerms := ionizations.MulPointwise(expBrace)

		fullMIntegral, err := minusTerms.TrapezoidIntegration(m.Parameters.CathodeFallLength, m.Parameters.GapLength)
		if err != nil {
			return nil, nil, err
		}
		exp2Lexp2dFullMIntegral := exp2L * exp2d * fullMIntegral

		fullPIntegral, err := plusTerms.TrapezoidIntegration(m.Parameters.CathodeFallLength, m.Parameters.GapLength)
		if err != nil {
			return nil, nil, err
		}

		for i := range density.Values { //for i := dIndex; i < len(density); i++ {
			x := m.XStep * float64(i)
			exp2x := expBrace(2 * x)
			dToxPIntegral, err := plusTerms.TrapezoidIntegration(m.Parameters.CathodeFallLength, x)
			if err != nil {
				return nil, nil, err
			}
			xToLPIntegral, err := plusTerms.TrapezoidIntegration(x, m.Parameters.GapLength)
			if err != nil {
				return nil, nil, err
			}
			dToxMIntegral, err := minusTerms.TrapezoidIntegration(m.Parameters.CathodeFallLength, x)
			if err != nil {
				return nil, nil, err
			}
			xToLMIntegral, err := minusTerms.TrapezoidIntegration(x, m.Parameters.GapLength)
			if err != nil {
				return nil, nil, err
			}
			density.Values[i] = max(0, scaleFactor*expBrace(-x)/normDivisor*utils.SumFloat64Slice([]float64{
				exp2Lexp2dFullMIntegral,
				exp2x * fullPIntegral,
				-exp2L * dToxPIntegral, //(utils.SumFloat64Slice(plusTerms[:i]) - utils.SumFloat64Slice(plusTerms[:dIndex])),
				-exp2d * xToLPIntegral,
				-exp2x * exp2d * dToxMIntegral,
				-exp2x * exp2L * xToLMIntegral,
			}))
			if math.IsInf(density.Values[i], 0) {
				panic("n(x) is Inf")
			}
		}
	}
	return []string{DensityPerCathodeElectronFluxKey}, []any{density}, nil
}

// ambipolar diffusion equation being solved numerically
func GlowDischargeDensity(m *model.Model) ([]string, []any, error) {
	densityPerCathodeElectronFlux := m.GetMetrics(model.DataHubKeyType(DensityPerCathodeElectronFluxKey)).(utils.GriddedInterval)
	firstDensityPeak := utils.ArgFirstPeak(densityPerCathodeElectronFlux.Values)
	densityPerCathodeElectronFlux.Values = append(densityPerCathodeElectronFlux.Values, densityPerCathodeElectronFlux.Values[len(densityPerCathodeElectronFlux.Values)-1])

	sumPerElectron := make([]int, m.TotalElectronsPassed)
	for e := range sumPerElectron {
		for x := range firstDensityPeak {
			sumPerElectron[e] += int(m.CollisionAtCell[lxgata.IONIZATION][x][e])
		}
		// j/q = phi_{e0} * (1 + \int_0^x_m { S_n(x) } dx)
		sumPerElectron[e] += 1
	}
	sumPlus1Mean, sumVariance := utils.MeanAndVariance(sumPerElectron, true)

	diffusionLoss := 0.
	if !m.Parameters.DoNotAccountDiffusionLoss {
		xm := densityPerCathodeElectronFlux.Step * float64(firstDensityPeak)
		d2MDensityIntegral, err := densityPerCathodeElectronFlux.TrapezoidIntegration(m.Parameters.CathodeFallLength, xm)
		if err != nil {
			return nil, nil, err
		}
		Lambda := m.Parameters.CathodeRadius / 2.4
		ambipolarDiffusionCoefficient := m.Parameters.IonMobility * m.Parameters.SlowElectronTemperature
		diffusionLoss = ambipolarDiffusionCoefficient * d2MDensityIntegral / (Lambda * Lambda)

	}
	sumPlus1Mean -= diffusionLoss

	//now sumPlus1Mean is 1 + 1/gamma with respect to abmipolar diffusion losses, that is 1 + integral(S(x) - D_a*n(x)/Lambda^2, 0, x0)

	sumMargin := utils.StudentedMargin(0.95, sumVariance, m.TotalElectronsPassed)
	cathodeElectronFlux := m.Parameters.CathodeCurrentDensity / (constants.ElementaryCharge * sumPlus1Mean)
	cathodeElectronFluxMargin := m.Parameters.CathodeCurrentDensity * sumMargin / (constants.ElementaryCharge * sumPlus1Mean * sumPlus1Mean)
	density := densityPerCathodeElectronFlux.Scale(cathodeElectronFlux)

	for e := range sumPerElectron {
		sumPerElectron[e] -= 1 // probably unnecessary measure against additive numeric error
	} //variance does not change after that act
	sumMean := utils.Mean(sumPerElectron) - diffusionLoss

	cathodeIonFlux := -cathodeElectronFlux * sumMean
	cathodeIonFluxMargin := cathodeElectronFluxMargin * sumMean

	return []string{
			GlowDischargeDensityKey,
			CathodeElectronFluxKey,
			CathodeElectronFluxMarginKey,
			CathodeIonFluxKey,
			CathodeIonFluxMarginKey,
			FirstDensityPeakIndexKey,
			SourceIntegralPerCathodeElectronFluxKey,
			SourceIntegralPerCathodeElectronFluxVarianceKey,
		}, []any{
			density,
			cathodeElectronFlux,
			cathodeElectronFluxMargin,
			cathodeIonFlux,
			cathodeIonFluxMargin,
			firstDensityPeak,
			sumMean,
			sumVariance,
		}, nil
}
