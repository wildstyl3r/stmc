package model

import "math"

// g
func (s *Model) VfromL(l float64) (V float64) {
	if l < s.Parameters.CathodeFallLength {
		//V(x) =  (-((V_c/d + C_E)(x/d - 2) + C_E)x + C_v)
		return math.FMA(
			-math.FMA(
				math.FMA(
					s.Vc,
					s.inverseCathodeFallLength,
					s.Parameters.ConstEField),
				math.FMA(
					l,
					s.inverseCathodeFallLength,
					-2.),
				s.Parameters.ConstEField),
			l,
			-(s.Vc + s.Va),
		)
	} else {
		// V(x) = -s.Va / (L-d) * (L - x) = C_E*(L-x)
		return s.Parameters.ConstEField * (s.Parameters.GapLength - l)
	}
}

// g^{-1}
func (s *Model) LfromV(V float64) (l float64) {
	if V < -s.Va {
		V_n := s.Parameters.ConstEField * s.Parameters.CathodeFallLength
		A := s.Vc + V_n
		l = s.Parameters.CathodeFallLength * (1 - (V_n+math.Sqrt(V_n*V_n-4*A*(V+s.Va)))/(2*A))
	} else {
		// V/C_E = L - x => x = L - V/C_E  ; [[---]=>[+-]=>[-]]
		return math.FMA(-V, s.inverseConstEField, s.Parameters.GapLength)
	}
	if math.IsNaN(l) {
		panic("x is NaN")
	}
	return
}

func (s *Model) EFieldFromL(l float64) (E float64) { // V/m
	if l < s.Parameters.CathodeFallLength {
		return math.FMA(2*math.FMA(s.Vc, s.inverseCathodeFallLength, s.Parameters.ConstEField), math.FMA(l, s.inverseCathodeFallLength, -1), s.Parameters.ConstEField)
	} else {
		return s.Parameters.ConstEField
	}
}

func (s *Model) EFieldFromPotential(V float64) (E float64) {
	if V < -s.Va {
		V_n := s.Parameters.ConstEField * s.Parameters.CathodeFallLength
		A := s.Vc + V_n
		return -math.Sqrt(V_n*V_n-4*A*(V+s.Va)) * s.inverseCathodeFallLength
	} else {
		return s.Parameters.ConstEField
	}

	// E = s.Parameters.ConstEField
	// if V < -s.Va {
	// 	E = math.FMA(-2*s.inverseCathodeFallLength, math.Sqrt(-(V+s.Va)*s.Vc), E)
	// }
	// return
}
