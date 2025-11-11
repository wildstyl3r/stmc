package utils

func PolynomialDerivative(x float64, coefficients []float64) (d float64) {
	derivCoefficients := make([]float64, len(coefficients)-1)
	for i := range coefficients {
		if i > 0 {
			derivCoefficients[i-1] = coefficients[i] * float64(i)
		}
	}
	for i := range derivCoefficients {
		d *= x
		d += derivCoefficients[len(derivCoefficients)-1-i]
	}
	return d
}
