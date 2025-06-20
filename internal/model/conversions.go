package model

import "math"

// g
func (s *Model) VfromL(l float64) (V float64) {
	cathodeFallPortion := l * s.inverseCathodeFallLength
	if l < s.Parameters.CathodeFallLength {
		V = s.Vc * (cathodeFallPortion*(2.-cathodeFallPortion) - 1.)
	} else {
		V = -s.Parameters.ConstEField * (l - s.Parameters.CathodeFallLength)
	}
	V -= s.Va
	return V
}

// g^{-1}
func (s *Model) LfromV(V float64) (l float64) {
	if V < -s.Va {
		l = math.FMA(-s.Parameters.CathodeFallLength, math.Sqrt((V+s.Va)*s.inverseNegativeVc), s.Parameters.CathodeFallLength)
	} else {
		l = math.FMA((V + s.Va), s.inverseAbsConstEField, s.Parameters.CathodeFallLength)
	}
	if math.IsNaN(l) {
		panic("x is NaN")
	}
	return
}

func (s *Model) EFieldFromL(l float64) (E float64) { // V/m
	E = s.Parameters.ConstEField
	if l < s.Parameters.CathodeFallLength {
		E = math.FMA(-2.*(1.-l*s.inverseCathodeFallLength)*s.Vc, s.inverseCathodeFallLength, E)
	}
	return
}

func (s *Model) EFieldFromPotential(V float64) (E float64) {
	E = s.Parameters.ConstEField
	if V < -s.Va {
		E = math.FMA(-2*s.inverseCathodeFallLength, math.Sqrt(-(V+s.Va)*s.Vc), E)
	}
	return
}
