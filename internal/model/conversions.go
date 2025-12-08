package model

import "math"

// g
func (s *Model) VfromL(x float64) (V float64) {
	if x < s.Parameters.CathodeFallLength { //&& r2 < s.cathodeRadius2 {
		//V(x) =  (-((V_c/d + C_E)(x/d - 2) + C_E)x + C_v)
		V = math.FMA(
			-math.FMA(
				math.FMA(
					s.Vc,
					s.inverseCathodeFallLength,
					s.Parameters.ConstEField),
				math.FMA(
					x,
					s.inverseCathodeFallLength,
					-2.),
				s.Parameters.ConstEField),
			x,
			-(s.Vc + s.Va),
		)
	} else { //if r2 > s.boundaryRadius2 {
		// V(x) = -s.Va / (L-d) * (L - x) = C_E*(L-x)
		V = s.Parameters.ConstEField * (s.Parameters.SimulationLength() - x)
		// } else {
		// 	return s.VfromL(x, 0) * (s.boundaryThickness + s.Parameters.CathodeRadius - math.Sqrt(r2)) / s.boundaryThickness
	}
	if math.IsNaN(V) {
		println("V is nan")
	}
	return
}

// g^{-1}
func (s *Model) LfromV(V float64) (l float64) {
	if V < -s.Va { //&& r2 < s.cathodeRadius2 {
		V_n := s.Parameters.ConstEField * s.Parameters.CathodeFallLength
		A := s.Vc + V_n
		l = s.Parameters.CathodeFallLength * (1 - (V_n+math.Sqrt(V_n*V_n-4*A*(V+s.Va)))/(2*A))
		// } else if r2 < s.boundaryRadius2 {
		// 	return s.LfromV(V*s.boundaryThickness/(s.boundaryThickness+s.Parameters.CathodeRadius-math.Sqrt(r2)), 0)
	} else {
		// V/C_E = L - x => x = L - V/C_E  ; [[---]=>[+-]=>[-]]
		return math.FMA(-V, s.inverseConstEField, s.Parameters.SimulationLength())
	}
	if math.IsNaN(l) {
		panic("x is NaN")
	}
	return
}

func (s *Model) EFieldFromL(x float64) (E float64) { // V/m
	if x < s.Parameters.CathodeFallLength { // && r2 < s.cathodeRadius2 {
		return math.FMA(2*math.FMA(s.Vc, s.inverseCathodeFallLength, s.Parameters.ConstEField), math.FMA(x, s.inverseCathodeFallLength, -1), s.Parameters.ConstEField)
		// } else if r2 < s.boundaryRadius2 {
		// 	radial := s.VfromL(x, 0) / s.boundaryThickness
		// 	return -math.Sqrt(radial*radial + s.Parameters.ConstEField*s.Parameters.ConstEField)
	} else {
		return s.Parameters.ConstEField
	}
}

func (s *Model) EFieldFromPotential(V float64) (E float64) {
	if V < -s.Va { //&& r2 < discharge {
		V_n := s.Parameters.ConstEField * s.Parameters.CathodeFallLength
		A := s.Vc + V_n
		return -math.Sqrt(V_n*V_n-4*A*(V+s.Va)) * s.inverseCathodeFallLength
		// } else if r2 < boundary {
		// 	radial := V / (s.boundaryThickness + s.Parameters.CathodeRadius - math.Sqrt(r2))
		// 	return -math.Sqrt(radial*radial + s.Parameters.ConstEField*s.Parameters.ConstEField)
	} else {
		return s.Parameters.ConstEField
	}
}

// func TimeToIntersectCircle(y0,z0,R,absV,eta float64) (t float64) {
// k,l,w := y0/absV,z0/absV,R/absV
// u := k*math.Cos(eta)+l*math.Sin(eta)
// D := 4*(u*u+w*w)
// }

// func (s *Model) PositionFromTime(p Particle, e0, t float64) (x float64) {

// }

// func (s *Model) EnergyFromTime(p Particle, e0, t float64) (e1 float64) {

// }
