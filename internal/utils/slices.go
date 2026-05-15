package utils

import (
	"cmp"
	"math"
	"slices"
	"sort"

	"golang.org/x/exp/constraints"
)

type MutableMap[K comparable, V any] struct {
	keys   map[K]int
	values []V
}

func NewMutableMap[K comparable, V any]() MutableMap[K, V] {
	return MutableMap[K, V]{
		keys:   make(map[K]int),
		values: make([]V, 0),
	}
}

func (m *MutableMap[K, V]) Keys() (keys []K) {
	for k := range m.keys {
		keys = append(keys, k)
	}
	return
}

func (m *MutableMap[K, V]) GetPointer(key K) *V {
	if index, exist := m.keys[key]; exist {
		return &m.values[index]
	} else {
		return nil
	}
}

func Argmax[T cmp.Ordered](arr []T) (argmax int) {
	for i := range arr {
		if cmp.Compare(arr[i], arr[argmax]) == 1 {
			argmax = i
		}
	}
	return
}

func ArgFirstPeak[T cmp.Ordered](arr []T) (argFirstPeak int) {
	for i := range arr {
		if i+2 < len(arr) && cmp.Compare(arr[i], arr[i+2]) == 1 && cmp.Compare(arr[i], arr[i+1]) == 1 && cmp.Compare(arr[i+1], arr[i+2]) == 1 {
			return i
		}
	}
	return len(arr) - 1
}

func Max[T cmp.Ordered](arr []T) (max T) {
	return arr[Argmax(arr)]
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

func ArgsortIndexable(s []ResultInterface, abs bool) (indices []int) {
	indices = make([]int, len(s))
	for i := range indices {
		indices[i] = i
	}

	sort.Slice(indices, func(i, j int) bool {
		return cmp.Compare(s[indices[i]].Index(), s[indices[j]].Index()) == -1
	})

	return indices
}

func SumIntSlice[T constraints.Integer](arr []T) (sum T) {
	for i := range arr {
		sum += arr[i]
	}
	return
}

func SumFloat64Slice(arr []float64, presort bool) (sum float64) { // Kahan's algorithm
	compensation := 0.
	if presort {
		summationOrder := Argsort(arr, true)
		for i := range summationOrder {
			y := arr[summationOrder[i]] - compensation
			temp := sum + y
			compensation = (temp - sum) - y
			sum = temp
		}
	} else {
		for i := range arr {
			y := arr[i] - compensation
			temp := sum + y
			compensation = (temp - sum) - y
			sum = temp
		}
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

func AnyIntersection(a, b []string) string {
	for i := range a {
		if slices.Contains(b, a[i]) {
			return a[i]
		}
	}
	return ""
}

func FullIntersection(a, b []string) (intersection []string) {
	for i := range a {
		if slices.Contains(b, a[i]) {
			intersection = append(intersection, a[i])
		}
	}
	return intersection
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

func NonInfIndicies(s []float64) (ixs []int) {
	for i, v := range s {
		if !math.IsInf(v, 0) {
			ixs = append(ixs, i)
		}
	}
	return
}

func Select[T any](s []T, ixs []int) (t []T) {
	for _, i := range ixs {
		t = append(t, s[i])
	}
	return
}

func ExcludeView[T any](i int, s []T) []T {
	if i >= len(s) {
		return s
	}
	newSlice := make([]T, len(s)-1)
	copy(newSlice, s[:i])
	copy(newSlice[i:], s[i+1:])
	return newSlice
}
