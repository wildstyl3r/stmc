package model

import (
	"math"
)

// g
func (s *Model) VfromL(x float64) (V float64) {
	if !s.Parameters.DarkDischarge && x < s.FieldSwitchPosition { //&& r2 < s.cathodeRadius2 {
		//V(x) =  2V/d * x * (1 - x/(2d)) + V_OFFSET = [-2V/d * x | * ([x | /(2d) | - 1]) | + V_OFFSET ]
		//V = math.FMA(s.EFieldB*x, math.FMA(x, 0.5*s.inverseCathodeFallLength, -1), s.VOffset)
		V = -s.Vs*math.Pow(math.FMA(-x, s.inverseCathodeFallLength, 1), s.Parameters.GetSheathFieldPower()+1.) + s.SheathPotentialCorrection
	} else {
		// V(x) = -s.Va / (L-d) * (L - x) = C_E*(L-x)
		V = s.Parameters.ConstEField * (s.Parameters.SimulationLength() - x)
	}
	if math.IsNaN(V) {
		println("V is nan")
	}
	return
}

// g^{-1}
func (s *Model) LfromV(V float64) (l float64) {
	if !s.Parameters.DarkDischarge && V < s.FieldSwitchPotential {
		l = s.Parameters.CathodeFallLength * (1 - math.Pow((s.SheathPotentialCorrection-V)*s.inverseCathodePotential, 1./(s.Parameters.GetSheathFieldPower()+1)))
		//l = s.Parameters.CathodeFallLength * (1 - math.Sqrt((s.VLSkip-V)/s.Vs))
	} else {
		// V/C_E = L - x => x = L - V/C_E  ; [[---]=>[+-]=>[-]]
		return math.FMA(-V, s.inverseConstEField, s.Parameters.SimulationLength())
	}
	if math.IsNaN(l) {
		panic("x is NaN")
	}
	return
}

func (s *Model) LNodefromVCached(V float64) (i int) {
	if !s.Parameters.DarkDischarge && V < s.FieldSwitchPotential-5 {
		i = s.lookupCellNodeFromPotential[int(-V/s.Parameters.EnergyStep)]
		//l = s.Parameters.CathodeFallLength * (1 - math.Pow((s.SheathPotentialCorrection-V)*s.inverseCathodePotential, 1./(s.Parameters.GetSheathFieldPower()+1)))
		//l = s.Parameters.CathodeFallLength * (1 - math.Sqrt((s.VLSkip-V)/s.Vs))
	} else if !s.Parameters.DarkDischarge && V < s.FieldSwitchPotential {
		i = int(s.Parameters.CathodeFallLength * (1 - math.Pow((s.SheathPotentialCorrection-V)*s.inverseCathodePotential, 1./(s.Parameters.GetSheathFieldPower()+1))) / s.XStep)
	} else {
		// V/C_E = L - x => x = L - V/C_E  ; [[---]=>[+-]=>[-]]
		return int(math.FMA(-V, s.inverseConstEField, s.Parameters.SimulationLength()) / s.XStep)
	}
	return
}

func (s *Model) EFieldFromL(x float64) (E float64) { // V/m
	if !s.Parameters.DarkDischarge && x < s.FieldSwitchPosition { // && r2 < s.cathodeRadius2 {
		return s.EFieldScale * math.Pow(math.FMA(-x, s.inverseCathodeFallLength, 1), s.Parameters.GetSheathFieldPower())
		// } else if r2 < s.boundaryRadius2 {
		// 	radial := s.VfromL(x, 0) / s.boundaryThickness
		// 	return -math.Sqrt(radial*radial + s.Parameters.ConstEField*s.Parameters.ConstEField)
	} else {
		return s.Parameters.ConstEField
	}
}

func (s *Model) EFieldFromPotential(V float64) (E float64) {
	if !s.Parameters.DarkDischarge && V < s.FieldSwitchPotential { //&& r2 < discharge {
		k := s.Parameters.GetSheathFieldPower()
		return s.EFieldScale * math.Pow((s.SheathPotentialCorrection-V)*s.inverseCathodePotential, k/(k+1))
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
