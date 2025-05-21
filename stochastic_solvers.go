package main

import (
	"fmt"
	"math/rand"
)

func stochasticApproximation(leftBound, rightBound,
	initialTheta,
	approxFDerivative,
	thetaPrecision float64, averagingLastN int,
	f func(float64) float64) float64 {
	thetas := []float64{initialTheta}
	fValues := []float64{f(thetas[0])}
	var a = 1 / approxFDerivative //((fRight - fLeft) / (right - left))
	fmt.Printf("calculated a-factor value: %v\n", a)
	for i := 0; i < averagingLastN || diameterOfSet(thetas[max(i-averagingLastN, 0):]) > thetaPrecision; i++ {
		newTheta := thetas[i] - fValues[i]*a/float64(i+1)
		if newTheta < leftBound || rightBound < newTheta {
			newTheta = leftBound + rand.Float64()*(rightBound-leftBound)
			println("RAND theta")
		}
		thetas = append(thetas, newTheta)
		fValues = append(fValues, f(newTheta))
	}
	return thetas[len(thetas)-1]
}
