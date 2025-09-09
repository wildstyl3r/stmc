package utils

import (
	"cmp"
	"math"
	"slices"
	"sort"

	"golang.org/x/exp/constraints"
)

func Argmax[T cmp.Ordered](arr []T) (argmax int) {
	for i := range arr {
		if cmp.Compare(arr[i], arr[argmax]) == 1 {
			argmax = i
		}
	}
	return
}

func Argsort(s []float64, abs bool) (indices []int) {
	indices = make([]int, len(s))
	for i := range indices {
		indices[i] = i
	}

	if abs {
		sort.Slice(indices, func(i, j int) bool {
			return math.Abs(s[indices[i]]) < math.Abs(s[indices[j]])
		})
	} else {
		sort.Slice(indices, func(i, j int) bool {
			return s[indices[i]] < s[indices[j]]
		})
	}

	return indices
}

func SumIntSlice[T constraints.Integer](arr []T) (sum T) {
	for i := range arr {
		sum += arr[i]
	}
	return
}

func SumFloat64Slice(arr []float64) (sum float64) { // Kahan's algorithm
	compensation := 0.
	summationOrder := Argsort(arr, true)
	for i := range summationOrder {
		y := arr[summationOrder[i]] - compensation
		temp := sum + y
		compensation = (temp - sum) - y
		sum = temp
	}
	return sum
}

func Apply[S, T any](f func(S) T, data []S) (result []T) {
	result = make([]T, len(data))
	for i := range result {
		result[i] = f(data[i])
	}
	return
}

func Intersect(a, b []string) *string {
	for i := range a {
		if slices.Contains(b, a[i]) {
			return &a[i]
		}
	}
	return nil
}

type constErr string

const tooManySubsets constErr = "the set is too large (> 16 elements) for subset enumeration"

func (ce constErr) Error() string {
	return string(ce)
}

func PartitionByMask[T any](s []T, mask uint16) (first, second []T, err error) {
	if len(s) > 16 {
		return nil, nil, tooManySubsets
	}
	if len(s) == 0 {
		return nil, nil, nil
	}
	for j := range s {
		if mask&(1<<j) == 0 {
			first = append(first, s[j])
		} else {
			second = append(second, s[j])
		}
	}
	return
}

func DiameterOfSet(s []float64) (d float64) {
	for i := range s {
		for j := i + 1; j < len(s); j++ {
			d = max(d, math.Abs(s[i]-s[j]))
		}
	}
	return
}
