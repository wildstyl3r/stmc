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
				fMargin = StudentedMargin(confidence, fValues[max(len(fValues)-minSteps, len(fValues)/4):])
				meanF = Mean(fValues[max(len(fValues)-minSteps, len(fValues)/4):])
				fmt.Printf("SA's mean function value is %f+/-%f required to be around 0 with margin less than %f \n", meanF, fMargin, precision)
			}
		} else {
			if len(thetas) > minSteps {
				thetaMargin = StudentedMargin(confidence, thetas[max(len(thetas)-minSteps, len(thetas)/4):])
				println("SA's theta confidence interval is ", thetaMargin*2, " required to be less than ", precision)
				// fConfidenceInterval = StudentedConfidenceInterval(confidence, fValues[len(thetas)/3:])
			}
		}
	}
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
			thetaMargin = Vec2D{StudentedMargin(confidence, xs), StudentedMargin(confidence, ys)}
			fmt.Printf("SA's theta confidence interval is %v required to be less than %v", thetaMargin, thetaPrecision)
		}
	}
	xs, ys := Vecs2Dto2slices(thetas[len(thetas)-minSteps:])
	return Vec2D{Mean(xs), Mean(ys)}, thetaMargin
}
