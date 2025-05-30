package sec

import (
	"fmt"
	"math"
	"math/rand"
)

func stochasticApproximation(leftBound, rightBound,
	initialTheta,
	approxFDerivative,
	thetaPrecision,
	quantile float64, minSteps int,
	f func(float64) float64) float64 {
	thetas := []float64{initialTheta}
	sumThetas := initialTheta
	sumSquareThetas := initialTheta * initialTheta
	fValues := []float64{f(thetas[0])}
	var a = 1 / approxFDerivative //((fRight - fLeft) / (right - left))
	fmt.Printf("\ncalculated a-factor value: %v\n\n", a)
	var confidenceInterval float64 = math.Inf(1)
	for i := 0; confidenceInterval > thetaPrecision; i++ {
		newTheta := thetas[i] - fValues[i]*a/float64(i+1)
		if newTheta < leftBound || rightBound < newTheta {
			newTheta = leftBound + rand.Float64()*(rightBound-leftBound)
			println("RAND theta")
		}
		thetas = append(thetas, newTheta)
		fValues = append(fValues, f(newTheta))

		sumThetas += newTheta
		sumSquareThetas += newTheta * newTheta

		if i >= minSteps {
			sumThetas -= thetas[i-minSteps]
			sumSquareThetas -= thetas[i-minSteps] * thetas[i-minSteps]
		}

		if len(thetas) > minSteps {
			varianceThetaNow := variance(thetas[len(thetas)-minSteps:], true)
			confidenceInterval = 2 * math.Sqrt(varianceThetaNow) * quantile / math.Sqrt(float64(minSteps))
		}
	}
	return sumThetas / float64(len(thetas))
}
