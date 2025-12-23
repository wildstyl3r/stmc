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
	tests := []struct{ x, e float64 }{
		{0, 0},
		{0.0002, 15},
		{0.00032, 5},
		{0.000032, 245},
	}
	for i := range tests {
		particle := Particle{
			x:        tests[i].x,
			eKinetic: tests[i].e,
		}
		particle.recalcParams(&m)
		checkX, _ := particle.trajectory.PositionAfterTime(particle.x, 0, &m, true)
		if math.Abs(checkX-tests[i].x) > 1e-9 {
			t.Errorf("want x0 = %f, got %f", tests[i].x, checkX)
			return
		}
		checkX, _ = particle.trajectory.PositionAfterTime(particle.x, particle.trajectory.TimeBetweenPositionsNoReversal(particle.x, 4*particle.x, &m, true), &m, true)
		if math.Abs(checkX-particle.x*4) > 1e-9 {
			t.Errorf("want x0 = %f, got %f", particle.x*4, checkX)
			return
		}
	}

}
