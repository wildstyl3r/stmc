package model

import (
	"math"
	"testing"

	"github.com/wildstyl3r/stmc/internal/config"
)

func TestEFromL(t *testing.T) {
	tests := []*Model{
		{
			Parameters: config.ModelParameters{
				CathodeFallLength:    0.001,
				CathodeFallPotential: 1000,
				ConstEField:          -100,
				GapLength:            0.005,
			}},
	}
	for i := range tests {
		tests[i].InitVars()
		p := tests[i]
		step := p.Parameters.CathodeFallLength / 200.
		for xIndex := range 1000 {
			x := float64(xIndex) * step
			convE := p.EFieldFromL(x)
			var fullE float64
			if x < p.Parameters.CathodeFallLength {
				fullE = -2*p.Parameters.CathodeFallPotential/p.Parameters.CathodeFallLength*(1-x/p.Parameters.CathodeFallLength) + p.Parameters.ConstEField
			} else {
				fullE = p.Parameters.ConstEField
			}

			if math.Abs(convE-fullE) > 1e-5 {
				t.Errorf("want E(x=%f) = %f, got %f", x, fullE, convE)
			}
		}
		cathodeConvV := p.VfromL(0)
		fullConvV := -p.Parameters.CathodeFallPotential + p.Parameters.ConstEField*p.Parameters.SimulationLength()
		if math.Abs(cathodeConvV-fullConvV) > 1e-5 {
			t.Errorf("want V(x=0) = %f, got %f", fullConvV, cathodeConvV)
		}
		cathodeConvV = p.VfromL(0.5 * p.Parameters.CathodeFallLength)
		fullConvV = -(1-0.75)*p.Parameters.CathodeFallPotential + p.Parameters.ConstEField*(p.Parameters.SimulationLength()-0.5*p.Parameters.CathodeFallLength)
		if math.Abs(cathodeConvV-fullConvV) > 1e-5 {
			t.Errorf("want V(x=%f) = %f, got %f", p.Parameters.CathodeFallLength*0.5, fullConvV, cathodeConvV)
		}
		cathodeConvV = p.VfromL(0.2 * p.Parameters.CathodeFallLength)
		fullConvV = -(1-0.36)*p.Parameters.CathodeFallPotential + p.Parameters.ConstEField*(p.Parameters.SimulationLength()-0.2*p.Parameters.CathodeFallLength)
		if math.Abs(cathodeConvV-fullConvV) > 1e-5 {
			t.Errorf("want V(x=%f) = %f, got %f", p.Parameters.CathodeFallLength*0.2, fullConvV, cathodeConvV)
		}
		cathodeConvV = p.VfromL(p.Parameters.CathodeFallLength)
		fullConvV = p.Parameters.ConstEField * (p.Parameters.SimulationLength() - p.Parameters.CathodeFallLength)
		if math.Abs(cathodeConvV-fullConvV) > 1e-5 {
			t.Errorf("want V(x=%f) = %f, got %f", p.Parameters.CathodeFallLength, fullConvV, cathodeConvV)
		}
	}
}
