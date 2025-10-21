package utils

import (
	"math"
	"testing"
)

func TestVariableLimitedDoubleIntegral(t *testing.T) {
	nodes := 100
	quarter := make([]float64, nodes+1)
	for i := range quarter {
		quarter[i] = 4
	}
	step := 1. / float64(nodes)
	interval := GriddedInterval{Values: quarter, Step: step}
	value, err := interval.VariableLimitDoubleIntegration(0, 1, step, 0, func(x float64) float64 { return math.Sqrt(1 - x*x) })
	t.Logf("number of nodes is %d, value is %f, error is %v", nodes, value, err)
	if err != nil {
		t.Errorf("error: %v", err)
		return
	}
	if math.Abs(value-math.Pi) > step {
		t.Errorf("expected %f, got %f", math.Pi, value)
	}
}
