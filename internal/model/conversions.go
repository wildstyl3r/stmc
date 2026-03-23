package model

import (
	"math"

	"github.com/wildstyl3r/stmc/internal/constants"
	"github.com/wildstyl3r/stmc/internal/utils"
)

// g
func (s *Model) VfromL(x float64) (V float64) {
	if x < s.FieldSwitchPosition { //&& r2 < s.cathodeRadius2 {
		//V(x) =  2V/d * x * (1 - x/(2d)) + V_OFFSET = [-2V/d * x | * ([x | /(2d) | - 1]) | + V_OFFSET ]
		V = math.FMA(s.EFieldB*x, math.FMA(x, 0.5*s.inverseCathodeFallLength, -1), s.VOffset)
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
	if V < s.FieldSwitchPotential {
		l = s.Parameters.CathodeFallLength * (1 - math.Sqrt((s.VLSkip-V)/s.Vs))
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
	if x < s.FieldSwitchPosition { // && r2 < s.cathodeRadius2 {
		return math.FMA(s.EFieldA, x, s.EFieldB)
		// } else if r2 < s.boundaryRadius2 {
		// 	radial := s.VfromL(x, 0) / s.boundaryThickness
		// 	return -math.Sqrt(radial*radial + s.Parameters.ConstEField*s.Parameters.ConstEField)
	} else {
		return s.Parameters.ConstEField
	}
}

func (s *Model) EFieldFromPotential(V float64) (E float64) {
	if V < s.FieldSwitchPotential { //&& r2 < discharge {
		return s.EFieldB * math.Sqrt((s.VLSkip-V)/s.Vs)
		// } else if r2 < boundary {
		// 	radial := V / (s.boundaryThickness + s.Parameters.CathodeRadius - math.Sqrt(r2))
		// 	return -math.Sqrt(radial*radial + s.Parameters.ConstEField*s.Parameters.ConstEField)
	} else {
		return s.Parameters.ConstEField
	}
}

func (m *Model) RateIncrement(e1, e2, totalEnergy float64, ks []utils.KahanSummable) {
	if math.Abs(e1-e2) < 1e-9 {
		return
	}
	cs1 := m.IonizationCSLookup[int(e1/m.Parameters.EnergyStep)] / -constants.ElementaryCharge
	cs2 := m.IonizationCSLookup[int(e2/m.Parameters.EnergyStep)] / -constants.ElementaryCharge
	if cs2 == 0 {
		return
	}
	V1, V2 := totalEnergy-e1, totalEnergy-e2
	// ef1 := m.FinegrainedFieldFromPotentialLookup[int(V1/m.Parameters.EnergyStep/100)]
	// x1 := m.FinegrainedPositionFromPotentialLookup[int(V1/m.Parameters.EnergyStep/100)]
	// ef2 := m.FinegrainedFieldFromPotentialLookup[int(V2/m.Parameters.EnergyStep/100)]
	// x2 := m.FinegrainedPositionFromPotentialLookup[int(V2/m.Parameters.EnergyStep/100)]
	// if x2 == x1 {
	// 	return
	// }
	x1, x2 := m.LfromV(-V1), m.LfromV(-V2)
	ef1, ef2 := m.EFieldFromPotential(-V1), m.EFieldFromPotential(-V2)
	a, b := utils.GetLinearCoefficients(e1, cs1, e2, cs2)
	c, d := utils.GetLinearCoefficients(e1, ef1, e2, ef2)
	var increment float64
	if c == 0 {
		increment = (0.5*a/d*(e2+e1) + b/d) * (e2 - e1)
	} else {
		increment = ((b*c-a*d)/(c*c)*math.Log((c*e1+d)/(c*e2+d)) + a/c*(e2-e1))
	}
	dx := x2 - x1
	increment /= dx //m.Parameters.GasDensity
	x1Index, x2Index := int(math.Floor(x1/m.XStep)), int(math.Floor(x2/m.XStep))
	if x1Index >= 0 && x1Index < m.NumCells && x2Index >= 0 && x2Index < m.NumCells {
		for i := int(math.Ceil(x1 / m.XStep)); i < x2Index; i++ {
			m.TrajAveragedIonization[i].Add(increment)
		}
		m.TrajAveragedIonization[x1Index].Add((x1 - float64(x1Index)*m.XStep) / m.XStep * increment)
		m.TrajAveragedIonization[x2Index].Add((x2 - float64(x2Index)*m.XStep) / m.XStep * increment)
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
