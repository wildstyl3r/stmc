package model

import (
	"math"
	"testing"

	"github.com/wildstyl3r/stmc/internal/config"
)

func TestGetTimeMotionConstants(t *testing.T) {
	m := Model{
		Parameters: config.ModelParameters{
			CathodeFallLength:    0.001,
			CathodeFallPotential: 1000,
			ConstEField:          -100,
		},
		inverseCathodeFallLength: 1000,
	}
	m.CalculateStaticTimeMotionConstants()
	tests := []struct{ x, e float64 }{
		{0, 0},
		{0.0002, 15},
		{0.00032, 5},
		{0.000032, 245},
	}
	for i := range tests {
		c1, c2 := m.GetDynamicTimeMotionConstants(tests[i].x, tests[i].e, 1)
		checkX := m.SheathPositionAtTime(0, c1, c2)
		if math.Abs(checkX-tests[i].x) > 1e-9 {
			t.Errorf("want x0 = %f, got %f", tests[i].x, checkX)
			return
		}
		checkE := m.SheathEnergyAtWTime(0, c1, c2)
		if math.Abs(checkE-tests[i].e) > 1e-9 {
			t.Errorf("want e0 = %f, got %f", tests[i].e, checkE)
			return
		}
	}

}
