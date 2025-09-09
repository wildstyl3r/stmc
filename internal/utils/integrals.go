package utils

import (
	"fmt"
	"math"
)

type GriddedInterval struct {
	Step   float64
	Offset float64
	Values []float64
}

type IntervalError int

func (ie IntervalError) Error() string {
	return fmt.Sprintf("interval error: requested grid node %d is beyond interval", int(ie))
}

func (gi *GriddedInterval) getClosestNodeAbove(x float64) (node int, distance float64, err error) {
	x -= gi.Offset
	node = int(math.Ceil(x / gi.Step))
	if node < 0 {
		return node, math.Inf(-1), IntervalError(node)
	}
	if node >= len(gi.Values) {
		return len(gi.Values) - 1, x - float64(len(gi.Values))*gi.Step, IntervalError(node)
	}
	return node, float64(node)*gi.Step - x, nil
}

func (gi *GriddedInterval) getClosestNodeBelow(x float64) (node int, distance float64, err error) {
	x -= gi.Offset
	node = int(math.Floor(x / gi.Step))
	if node < 0 {
		return node, math.Inf(-1), IntervalError(node)
	}
	if node >= len(gi.Values) {
		return len(gi.Values) - 1, x - float64(len(gi.Values))*gi.Step, IntervalError(node)
	}
	return node, x - float64(node)*gi.Step, nil
}

func (gi *GriddedInterval) interpolate(x float64) (value float64, closestNodeBelow int, err error) {
	below, distance, err := gi.getClosestNodeBelow(x)
	if err != nil {
		return 0, below, err
	}
	exactNodeLocation := gi.Step * float64(below)
	if exactNodeLocation == x {
		return gi.Values[below], below, nil
	}
	above := below + 1
	if above >= len(gi.Values) {
		return 0, below, IntervalError(above)
	}
	valueBelow, valueAbove := gi.Values[below], gi.Values[above]
	return valueBelow + (valueAbove-valueBelow)*(distance/gi.Step), below, nil

}

func (gi *GriddedInterval) Interpolate(x float64) float64 {
	value, _, err := gi.interpolate(x)
	if err != nil {
		fmt.Printf("error interpolating %v at %f: %#v", gi, x, err)
	}
	return value
}

func TrapezoidIntegration(lowerLimit, upperLimit float64, f GriddedInterval) (integral float64, err error) {
	if lowerLimit == upperLimit {
		return 0, nil
	}
	negative := false
	if lowerLimit > upperLimit {
		negative = true
		lowerLimit, upperLimit = upperLimit, lowerLimit
	}

	valueAtLowerLimit, gridNodeBelowLowerLimitInt, err := f.interpolate(lowerLimit)
	if err != nil {
		return 0, err
	}
	gridNodeAboveLowerLimitInt := gridNodeBelowLowerLimitInt + 1
	valueAtUpperLimit, gridNodeBelowUpperLimitInt, err := f.interpolate(upperLimit)
	if err != nil {
		return 0, err
	}

	if gridNodeAboveLowerLimitInt > gridNodeBelowUpperLimitInt { // upper and lower limits are located inside the same cell
		deltaArg := upperLimit - lowerLimit
		integral = 0.5 * (valueAtLowerLimit + valueAtUpperLimit) * deltaArg
	} else {
		lowerNode, lowerStep, err := f.getClosestNodeAbove(lowerLimit)
		if err != nil {
			return 0, err
		}
		lowerPiece := 0.5 * (valueAtLowerLimit + f.Values[lowerNode]) * lowerStep

		upperNode, upperStep, err := f.getClosestNodeBelow(upperLimit)
		if err != nil {
			return 0, err
		}
		upperPiece := 0.5 * (valueAtUpperLimit + f.Values[upperNode]) * upperStep

		sumTerms := []float64{lowerPiece, upperPiece}
		for i := lowerNode; i < upperNode; i++ {
			sumTerms = append(sumTerms, 0.5*(f.Values[i]+f.Values[i+1])*f.Step)
		}
		integral = SumFloat64Slice(sumTerms)
	}

	if negative {
		integral = -integral
	}
	return
}

func VariableLimitDoubleIntegration(outerLowerLimit, outerUpperLimit, step, innerLowerLimit float64, innerUpperLimit func(float64) float64, f GriddedInterval) (integral float64, err error) { // addressing the case of int_a^b [ int_f(x)^g(x) dy ] dx
	negative := false
	if outerLowerLimit > outerUpperLimit {
		negative = true
		outerLowerLimit, outerUpperLimit = outerUpperLimit, outerLowerLimit
	}
	outerGrid := GriddedInterval{
		Values: make([]float64, int(math.Ceil((outerUpperLimit-outerLowerLimit)/step))+1),
		Step:   step,
		Offset: outerLowerLimit}
	for i := range outerGrid.Values {
		x := step*float64(i) + outerGrid.Offset
		outerGrid.Values[i], err = TrapezoidIntegration(innerLowerLimit, innerUpperLimit(x), f)
		if err != nil {
			return 0, err
		}
	}
	integral, err = TrapezoidIntegration(outerLowerLimit, outerUpperLimit, outerGrid)
	if negative {
		integral = -integral
	}
	return
}
