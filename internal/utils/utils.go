package utils

import (
	"cmp"
	"math"
	"math/rand"
	"slices"

	"golang.org/x/exp/constraints"

	"github.com/wildstyl3r/stmc/internal/constants"
)

func Argmax[T cmp.Ordered](arr []T) (argmax int) {
	for i := range arr {
		if cmp.Compare(arr[i], arr[argmax]) == 1 {
			argmax = i
		}
	}
	return
}

type Number interface {
	constraints.Float | constraints.Integer
}

func SumSlice[T Number](arr []T) (r T) {
	for i := range arr {
		r += arr[i]
	}
	return
}

func DiameterOfSet(s []float64) (d float64) {
	for i := range s {
		for j := i + 1; j < len(s); j++ {
			d = max(d, math.Abs(s[i]-s[j]))
		}
	}
	return
}

func Average[T Number](s []T) (mean float64) {
	for i := range s {
		mean += float64(s[i])
	}
	mean /= float64(len(s))
	return
}

func MeanAndVariance[T Number](s []T, unbiased bool) (mean, variance float64) {
	mean = Average(s)
	for i := range s {
		variance += (float64(s[i]) - mean) * (float64(s[i]) - mean)
	}
	if unbiased {
		variance /= float64(len(s) - 1)
	} else {
		variance /= float64(len(s))
	}

	return
}

func Variance[T Number](s []T, unbiased bool) float64 {
	_, v := MeanAndVariance(s, unbiased)
	return v
}

func IntAbs(a int) int {
	if a < 0 {
		return -a
	} else {
		return a
	}

}

func TableIntegrate(s []float64, multiply func(float64) float64, step float64) (sum float64) {
	for i := range s {
		if multiply == nil {
			sum += s[i]
		} else {
			sum += s[i] * multiply(float64(i)*step)
		}
	}
	sum *= step
	return
}

func R() float64 {
	return -math.Log(rand.Float64())
}

func UniformOnDisk(r float64) (a, b float64) {
	a, b = 2.*rand.Float64()-1., 2.*rand.Float64()-1.
	for a*a+b*b > 1. {
		a, b = 2.*rand.Float64()-1., 2.*rand.Float64()-1.
	}
	a *= r
	b *= r
	return
}

func Intersect(a, b []string) *string {
	for i := range a {
		if slices.Contains(b, a[i]) {
			return &a[i]
		}
	}
	return nil
}

func EV2J(val float64) float64 {
	return val * constants.ElectronCharge
}

func J2eV(val float64) float64 {
	return val / constants.ElectronCharge
}

func EV2electronVelocity(energy float64) (v float64) {
	v = math.Sqrt(2 * energy * constants.ElectronCharge / constants.ElectornMass)
	return
}
