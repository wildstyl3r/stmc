package extensions

import (
	"math"

	"github.com/wildstyl3r/lxgata"
	"github.com/wildstyl3r/stmc/internal/model"
	"github.com/wildstyl3r/stmc/internal/utils"
)

const (
	CathodeElectronCurrentFractionKey             = "CathodeElectronFlux"
	CathodeElectronCurrentFractionMarginKey       = "CathodeElectronFluxMargin"
	CathodeIonCurrentFractionKey                  = "CathodeIonFlux"
	CathodeIonCurrentFractionMarginKey            = "CathodeIonFluxMargin"
	FirstDensityPeakIndexKey                      = "FirstDensityPeakIndex"
	SourceIntegralPerCathodeElectronFluxKey       = "SingleElectronSourceIntegral"
	SourceIntegralPerCathodeElectronFluxMarginKey = "SingleElectronSourceIntegralMargin"
	DensityPerCathodeElectronFluxKey              = "DensityPerCathodeElectronFluxKey"
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

	if m.Parameters.NoDiffusionLoss {
		// n(x) = -int{int{S(h),d,chi},d, x} + c1*x + c2
		scaledIonizations := ionizations.Scale(1. / ambipolarDiffusionCoefficient)
		C1, err := scaledIonizations.RightTriangleIntegration(m.Parameters.CathodeFallLength, m.Parameters.GapLength, m.XStep, 0)
		if err != nil {
			return nil, nil, err
		}
		C1 /= (m.Parameters.GapLength - m.Parameters.CathodeFallLength)
		C2 := -C1 * m.Parameters.CathodeFallLength
		for i := range density.Values {
			x := float64(i) * m.XStep
			integral, err := scaledIonizations.RightTriangleIntegration(m.Parameters.CathodeFallLength, x, m.XStep, m.Parameters.CathodeFallLength)
			if err != nil {
				return nil, nil, err
			}
			density.Values[i] = -integral + C1*x + C2
		}
	} else {
		// n(x) = ( A / 2D ) * ([-x] / ([2d] - [2L]) ) * ( [2(L+d)]{d,L}m + [2x]{d,L}p - [2L]{d,x}p - [2d]{x,L}p - [2(x+d)]{d,x}m - [2(x+L)]{x,L}m )
		// n(x) = scale * [-x] / norm  * ( [2(L+d)]fullMinus + [2x]fullPlus - [2L]{d,x}p - [2d]{x,L}p - [2(x+d)]{d,x}m - [2(x+L)]{x,L}m )
		Lambda := m.GetMetrics(CharacteristicDiffusionScaleKey).(float64)
		// if m.Parameters.HyperbolicDensity {
		scaleFactor := Lambda / ambipolarDiffusionCoefficient

		constFactorBTerms := ionizations.MulPointwise(func(x float64) float64 { return math.Sinh((m.Parameters.GapLength - x) / Lambda) })
		// constFactorBTerms := utils.GriddedInterval{
		// 	Step:   m.XStep,
		// 	Offset: ionizations.Offset,
		// 	Values: make([]float64, m.NumCells)}
		// for i := range constFactorBTerms.Values {
		// 	x := m.XStep * float64(i)
		// 	constFactorBTerms.Values[i] = ionizations.Values[i] * math.Sinh((m.Parameters.GapLength-x)/Lambda)
		// }
		// constFactorBTerms.Values = append(constFactorBTerms.Values, constFactorBTerms.Values[len(constFactorBTerms.Values)-1])
		constFactorB, err := constFactorBTerms.TrapezoidIntegration(m.Parameters.CathodeFallLength, m.Parameters.GapLength, true)
		if err != nil {
			return nil, nil, err
		}
		constFactorB /= math.Sinh((m.Parameters.GapLength - m.Parameters.CathodeFallLength) / Lambda)

		convolution := utils.GriddedInterval{
			Step:   m.XStep,
			Offset: ionizations.Offset,
			Values: make([]float64, m.NumCells)}

		// P, Q, Pcompensation, Qcompensation := 0., 0., 0., 0.

		dIndex := int(m.Parameters.CathodeFallLength / m.XStep)

		for i := range density.Values {
			x := m.XStep * float64(i)
			if x > m.Parameters.CathodeFallLength {
				for j := dIndex; j <= i; j++ {
					xi := m.XStep * float64(j)
					convolution.Values[j] = ionizations.Values[j] * math.Sinh((x-xi)/Lambda)
				}
				conv, err := convolution.TrapezoidIntegration(m.Parameters.CathodeFallLength, x, false)
				if err != nil {
					return nil, nil, err
				}
				density.Values[i] = max(0, scaleFactor*utils.SumFloat64Slice([]float64{
					-conv,
					math.Sinh((x-m.Parameters.CathodeFallLength)/Lambda) * constFactorB,
				}, true))
			}
			// if x > m.Parameters.CathodeFallLength {
			// 	X := x / Lambda
			// 	{
			// 		y := ionizations.Values[i]*math.Sinh(X) - Pcompensation
			// 		temp := P + y
			// 		Pcompensation = (temp - P) - y
			// 		P = temp
			// 	}
			// 	{
			// 		y := ionizations.Values[i]*math.Cosh(X) - Qcompensation
			// 		temp := Q + y
			// 		Qcompensation = (temp - Q) - y
			// 		Q = temp
			// 	}
			// 	chX, shX := math.Cosh(X), math.Sinh(X)
			// 	density.Values[i] = max(0, scaleFactor*utils.SumFloat64Slice([]float64{
			// 		chX * P,
			// 		-shX * Q,
			// 		math.Sinh((x-m.Parameters.CathodeFallLength)/Lambda) * constFactorB,
			// 	}))
			// }

		}
		// } else {
		// 	powerFactor := 1. / Lambda
		// 	scaleFactor := Lambda / (2 * ambipolarDiffusionCoefficient)

		// 	expBrace := func(v float64) float64 { return math.Exp(v * powerFactor) }
		// 	exp2L, exp2d := expBrace(2*m.Parameters.GapLength), expBrace(2*m.Parameters.CathodeFallLength)
		// 	normDivisor := exp2d - exp2L

		// 	negExpBrace := func(x float64) float64 { return expBrace(-x) }
		// 	minusTerms := ionizations.MulPointwise(negExpBrace)
		// 	plusTerms := ionizations.MulPointwise(expBrace)

		// 	fullMIntegral, err := minusTerms.TrapezoidIntegration(m.Parameters.CathodeFallLength, m.Parameters.GapLength, true)
		// 	if err != nil {
		// 		return nil, nil, err
		// 	}
		// 	exp2Lexp2dFullMIntegral := exp2L * exp2d * fullMIntegral

		// 	fullPIntegral, err := plusTerms.TrapezoidIntegration(m.Parameters.CathodeFallLength, m.Parameters.GapLength, true)
		// 	if err != nil {
		// 		return nil, nil, err
		// 	}

		// 	for i := range density.Values { //for i := dIndex; i < len(density); i++ {
		// 		x := m.XStep * float64(i)
		// 		exp2x := expBrace(2 * x)
		// 		dToxPIntegral, err := plusTerms.TrapezoidIntegration(m.Parameters.CathodeFallLength, x, true)
		// 		if err != nil {
		// 			return nil, nil, err
		// 		}
		// 		xToLPIntegral, err := plusTerms.TrapezoidIntegration(x, m.Parameters.GapLength, true)
		// 		if err != nil {
		// 			return nil, nil, err
		// 		}
		// 		dToxMIntegral, err := minusTerms.TrapezoidIntegration(m.Parameters.CathodeFallLength, x, true)
		// 		if err != nil {
		// 			return nil, nil, err
		// 		}
		// 		xToLMIntegral, err := minusTerms.TrapezoidIntegration(x, m.Parameters.GapLength, true)
		// 		if err != nil {
		// 			return nil, nil, err
		// 		}
		// 		density.Values[i] = max(0, scaleFactor*expBrace(-x)/normDivisor*utils.SumFloat64Slice([]float64{
		// 			exp2Lexp2dFullMIntegral,
		// 			exp2x * fullPIntegral,
		// 			-exp2L * dToxPIntegral, //(utils.SumFloat64Slice(plusTerms[:i]) - utils.SumFloat64Slice(plusTerms[:dIndex])),
		// 			-exp2d * xToLPIntegral,
		// 			-exp2x * exp2d * dToxMIntegral,
		// 			-exp2x * exp2L * xToLMIntegral,
		// 		}, true))
		// 		if math.IsInf(density.Values[i], 0) {
		// 			panic("n(x) is Inf")
		// 		}
		// 		if math.IsNaN(density.Values[i]) {
		// 			println("n(x) is NaN")
		// 		}
		// 	}
		// }

	}
	return []string{DensityPerCathodeElectronFluxKey}, []any{density}, nil
}

// ambipolar diffusion equation being solved numerically
func GlowDischargeDensity(m *model.Model) ([]string, []any, error) {
	densityPerCathodeElectronFlux := m.GetMetrics(model.DataHubKeyType(DensityPerCathodeElectronFluxKey)).(utils.GriddedInterval)
	firstDensityPeak := utils.ArgFirstPeak(densityPerCathodeElectronFlux.Values)
	densityPerCathodeElectronFlux.Values = append(densityPerCathodeElectronFlux.Values, densityPerCathodeElectronFlux.Values[len(densityPerCathodeElectronFlux.Values)-1])
	// 	// j/q = phi_{e0} * (1 + \int_0^x_m { S_n(x) } dx)
	sumMean, sumMargin := m.IonizationsSumUpToCell[firstDensityPeak-1].MeanWithErrorMargin(0.95)
	aggPlus1 := m.IonizationsSumUpToCell[firstDensityPeak-1]
	aggPlus1.IncrementBy1()
	sumPlus1Mean := aggPlus1.Mean()

	diffusionLoss := 0.
	if !m.Parameters.NoDiffusionLoss {
		xm := densityPerCathodeElectronFlux.Step * float64(firstDensityPeak)
		d2MDensityIntegral, err := densityPerCathodeElectronFlux.TrapezoidIntegration(m.Parameters.CathodeFallLength, xm, true)
		if err != nil {
			return nil, nil, err
		}
		Lambda := m.GetMetrics(CharacteristicDiffusionScaleKey).(float64)
		ambipolarDiffusionCoefficient := m.GetMetrics(AmbipolarDiffusionCoefficientKey).(float64)
		diffusionLoss = ambipolarDiffusionCoefficient * d2MDensityIntegral / (Lambda * Lambda)

	}
	sumPlus1Mean -= diffusionLoss
	sumMean -= diffusionLoss

	//now sumPlus1Mean is 1 + 1/gamma with respect to abmipolar diffusion losses, that is 1 + integral(S(x) - D_a*n(x)/Lambda^2, 0, x0)
	var electronCurrentFraction, electronCurrentFractionMargin float64
	electronCurrentFraction = 1. / sumPlus1Mean
	electronCurrentFractionMargin = sumMargin / (sumPlus1Mean * sumPlus1Mean)

	// _sumMean, _sumVariance := m.IonizationsSumUpToCell[firstDensityPeak-1].MeanAndVariance()
	// fmt.Printf("aggregated sum mean and variance: [%f; %f], native: [%f, %f]\n", _sumMean, _sumVariance, utils.Mean(sumPerElectron), sumVariance)

	if math.IsNaN(sumMean) {
		println("ionization integral is nan")
	}

	cathodeIonCurrentFraction := sumMean / sumPlus1Mean
	cathodeIonCurrentFractionMargin := electronCurrentFractionMargin * sumMean

	return []string{
			CathodeElectronCurrentFractionKey,
			CathodeElectronCurrentFractionMarginKey,
			CathodeIonCurrentFractionKey,
			CathodeIonCurrentFractionMarginKey,
			FirstDensityPeakIndexKey,
			SourceIntegralPerCathodeElectronFluxKey,
			SourceIntegralPerCathodeElectronFluxMarginKey,
		}, []any{
			electronCurrentFraction,
			electronCurrentFractionMargin,
			cathodeIonCurrentFraction,
			cathodeIonCurrentFractionMargin,
			firstDensityPeak,
			sumMean,
			sumMargin,
		}, nil
}
