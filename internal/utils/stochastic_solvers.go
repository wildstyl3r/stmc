package utils

import (
	"fmt"
	"math"
	"math/rand"
)

func StochasticApproximation(lowerBound, upperBound,
	initialTheta,
	approxFDerivative,
	precision,
	// fPrecision,
	confidence float64, precisionForF bool, minSteps, maxSteps int,
	f func(float64) float64) (meanTheta, thetaMargin float64) {
	thetas := []float64{initialTheta}
	fValues := []float64{f(thetas[0])}
	var a = 1 / approxFDerivative //((fRight - fLeft) / (right - left))
	fmt.Printf("\ncalculated a-factor value: %v\n\n", a)
	// var fConfidenceInterval float64
	thetaMargin = math.Inf(1) //, fConfidenceInterval = , math.Inf(1) //&& fConfidenceInterval > fPrecision
	fMargin, meanF := math.Inf(1), math.Inf(1)
	for i := 0; ((!precisionForF && (thetaMargin*2 > precision)) || (precisionForF && (math.Abs(meanF)+fMargin > precision))) && i < maxSteps; i++ {
		newTheta := thetas[i] - fValues[i]*a/float64(i+1)
		if newTheta < lowerBound || upperBound < newTheta {
			newTheta = lowerBound + rand.Float64()*(upperBound-lowerBound)
			println("RAND theta")
		}
		thetas = append(thetas, newTheta)
		fValues = append(fValues, f(newTheta))

		if precisionForF {
			if len(fValues) > minSteps {
				fMargin = StudentedMarginFromData(confidence, fValues[max(len(fValues)-minSteps, len(fValues)/4):])
				meanF = Mean(fValues[max(len(fValues)-minSteps, len(fValues)/4):])
				fmt.Printf("SA's mean function value is %f+/-%f required to be around 0 with margin less than %f \n", meanF, fMargin, precision)
			}
		} else {
			if len(thetas) > minSteps {
				thetaMargin = StudentedMarginFromData(confidence, thetas[max(len(thetas)-minSteps, len(thetas)/4):])
				println("SA's theta confidence interval is ", thetaMargin*2, " required to be less than ", precision)
				// fConfidenceInterval = StudentedConfidenceInterval(confidence, fValues[len(thetas)/3:])
			}
		}
	}
	thetaMargin = StudentedMarginFromData(confidence, thetas[max(len(thetas)-minSteps, len(thetas)/4):])
	return Mean(thetas[max(len(thetas)-minSteps, len(thetas)/4):]), thetaMargin
}

func StochasticApproximationWithSameSizeSamples(lowerBound, upperBound,
	initialTheta,
	approxFDerivative,
	precision,
	// fPrecision,
	confidence float64, precisionForF bool, minSteps, maxSteps, sampleSize int,
	f func(float64) (meanEstimation float64, varianceEstimation float64), oracle func() *float64) (meanTheta, thetaMargin float64) {
	thetas := []float64{initialTheta}
	sMean0, sVar0 := f(thetas[0])
	fSMean, fSVar := []float64{sMean0}, []float64{sVar0}
	var a = 1 / approxFDerivative //((fRight - fLeft) / (right - left))
	fmt.Printf("\ncalculated a-factor value: %v\n\n", a)
	// var fConfidenceInterval float64
	thetaMargin = math.Inf(1) //, fConfidenceInterval = , math.Inf(1) //&& fConfidenceInterval > fPrecision
	fMargin, meanF := math.Inf(1), math.Inf(1)
	for i := 0; ((!precisionForF && (thetaMargin*2 > precision)) || (precisionForF && (math.Abs(meanF)+fMargin > precision))) && i < maxSteps; i++ {
		newTheta := thetas[i] - fSMean[i]/approxFDerivative/float64(i+1)
		suggestion := oracle()
		if suggestion != nil {
			newTheta = *suggestion
			fmt.Printf("accepting suggestion: %f\n", newTheta)
		}
		if newTheta < lowerBound || upperBound < newTheta {
			newTheta = lowerBound + rand.Float64()*(upperBound-lowerBound)
			println("RAND theta")
		}
		thetas = append(thetas, newTheta)
		sMean, sVar := f(newTheta)
		fSMean = append(fSMean, sMean)
		fSVar = append(fSVar, sVar)

		if i > 2 {
			trial := (0.5*approxFDerivative + 0.5*(fSMean[i-2]-fSMean[i-1])/(thetas[i-2]-thetas[i-1]))
			if !math.IsNaN(trial) {
				approxFDerivative = trial
			}
		}

		if precisionForF {
			if len(fSVar) > minSteps {
				subsampleLen := max(len(fSVar)-minSteps, len(fSVar)/3)
				fPoolVar := SameSizePoolVariance(fSVar[subsampleLen:])
				meanF = Mean(fSMean[subsampleLen:]) //SumFloat64Slice(fSMean[subsampleLen:]) / float64(subsampleLen)
				fMargin = StudentedMargin(confidence, fPoolVar, subsampleLen*sampleSize)
				fmt.Printf("SA's mean function value is %f+/-%f required to be around 0 with margin less than %f \n", meanF, fMargin, precision)
			}
		} else {
			if len(thetas) > minSteps {
				thetaMargin = StudentedMarginFromData(confidence, thetas[max(len(thetas)-minSteps, len(thetas)/3):])
				println("SA's theta confidence interval is ", thetaMargin*2, " required to be less than ", precision)
				// fConfidenceInterval = StudentedConfidenceInterval(confidence, fValues[len(thetas)/3:])
			}
		}
	}
	thetaMargin = StudentedMarginFromData(confidence, thetas[max(len(thetas)-minSteps, len(thetas)/4):])
	return Mean(thetas[max(len(thetas)-minSteps, len(thetas)/4):]), thetaMargin
}

func StochasticApproximation2D(lowerBound, upperBound,
	initialTheta Vec2D,
	approxJacobian Mat2D,
	thetaPrecision Vec2D,
	confidence float64, minSteps int,
	f func(Vec2D) Vec2D) (meanTheta, thetaMargin Vec2D) {
	thetas := []Vec2D{initialTheta}
	fValues := []Vec2D{f(thetas[0])}
	var a = approxJacobian.Inverse()
	fmt.Printf("\ncalculated a-factor value: %v\n\n", a)
	thetaMargin = Vec2D{math.Inf(1), math.Inf(1)}
	for i := 0; thetaMargin.X*2 > thetaPrecision.X && thetaMargin.Y*2 > thetaPrecision.Y; i++ {
		newTheta := thetas[i].Sub(a.MulVec(fValues[i]).Div(math.Sqrt(float64(i + 1))))
		if !newTheta.InsideOfRectangle(lowerBound, upperBound) {
			newTheta = lowerBound.Add(Vec2D{rand.Float64(), rand.Float64()}.PointMul(upperBound.Sub(lowerBound)))
			println("RAND theta")
		}
		thetas = append(thetas, newTheta)
		fValues = append(fValues, f(newTheta))

		if len(thetas) > minSteps {
			xs, ys := Vecs2Dto2slices(thetas[len(thetas)-minSteps:])
			thetaMargin = Vec2D{StudentedMarginFromData(confidence, xs), StudentedMarginFromData(confidence, ys)}
			fmt.Printf("SA's theta confidence interval is %v required to be less than %v", thetaMargin, thetaPrecision)
		}
	}
	xs, ys := Vecs2Dto2slices(thetas[len(thetas)-minSteps:])
	return Vec2D{Mean(xs), Mean(ys)}, thetaMargin
}
