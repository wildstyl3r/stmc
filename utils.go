package main

import (
	"bufio"
	"cmp"
	"encoding/csv"
	"fmt"
	"math"
	"math/rand"
	"os"
	"path/filepath"
	"slices"
	"sort"
	"strconv"
	"strings"

	"golang.org/x/exp/constraints"

	"github.com/facette/natsort"
	"github.com/wildstyl3r/lxgata"
)

const kBolzmann float64 = 1.380649e-23
const electronCharge = 1.602176634e-19                   // C
const electornMass float64 = 9.1093837139e-31            // [kg]
const Townsend float64 = 1.e-21                          // V * m^2
const freeSpacePermittivityE0 float64 = 8.8541878188e-12 // [m^-3 kg^{-1} s^4 A^2]
const quantile95 = 1.96

func ternarySearchMax(f func(float64) float64, left, right, eps float64) float64 {
	for right-left > eps {
		a := math.FMA(left, 2., right) / 3.
		b := math.FMA(right, 2., left) / 3.
		if f(a) > f(b) {
			right = b
		} else {
			left = a
		}
	}
	return (left + right) * 0.5
}

// return the point of the condition support that is not farther than eps from the support boundary
// invariant: at *right* condition must be TRUE
func binarySearch(condition func(float64) bool, falseDom, trueDom, eps float64) (float64, float64) {
	for math.Abs(trueDom-falseDom) > eps {
		c := (falseDom + trueDom) * 0.5
		if condition(c) {
			trueDom = c
		} else {
			falseDom = c
		}
	}
	return falseDom, trueDom
}

func ternarySearchMaxF(f func(float64) float64, left, right, eps float64) float64 {
	return f(ternarySearchMax(f, left, right, eps))
}

func argmax[T cmp.Ordered](arr []T) (argmax int) {
	for i := range arr {
		if cmp.Compare(arr[i], arr[argmax]) == 1 {
			argmax = i
		}
	}
	return
}

type Number interface {
	constraints.Float | constraints.Integer
}

func sumSlice[T Number](arr []T) (r T) {
	for i := range arr {
		r += arr[i]
	}
	return
}

func diameterOfSet(s []float64) (d float64) {
	for i := range s {
		for j := i + 1; j < len(s); j++ {
			d = max(d, math.Abs(s[i]-s[j]))
		}
	}
	return
}

func average[T Number](s []T) (mean float64) {
	for i := range s {
		mean += float64(s[i])
	}
	mean /= float64(len(s))
	return
}

func meanAndVariance[T Number](s []T, unbiased bool) (mean, variance float64) {
	mean = average(s)
	for i := range s {
		variance += (float64(s[i]) - mean) * (float64(s[i]) - mean)
	}
	if unbiased {
		variance /= float64(len(s) - 1)
	} else {
		variance /= float64(len(s))
	}

	return
}

func variance[T Number](s []T, unbiased bool) float64 {
	_, v := meanAndVariance(s, unbiased)
	return v
}

func intAbs(a int) int {
	if a < 0 {
		return -a
	} else {
		return a
	}

}

func tableIntegrate(s []float64, multiply func(float64) float64, step float64) (sum float64) {
	for i := range s {
		if multiply == nil {
			sum += s[i]
		} else {
			sum += s[i] * multiply(float64(i)*step)
		}
	}
	sum *= step
	return
}

func R() float64 {
	return -math.Log(rand.Float64())
}

func uniformOnDisk(r float64) (a, b float64) {
	a, b = 2.*rand.Float64()-1., 2.*rand.Float64()-1.
	for a*a+b*b > 1. {
		a, b = 2.*rand.Float64()-1., 2.*rand.Float64()-1.
	}
	a *= r
	b *= r
	return
}

func intersect(a, b []string) *string {
	for i := range a {
		if slices.Contains(b, a[i]) {
			return &a[i]
		}
	}
	return nil
}

func eV2J(val float64) float64 {
	return val * electronCharge
}

func J2eV(val float64) float64 {
	return val / electronCharge
}

func eV2electronVelocity(energy float64) (v float64) {
	v = math.Sqrt(2 * energy * electronCharge / electornMass)
	return
}

func readFloatPairs(filename string) ([][]float64, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, fmt.Errorf("error opening file: %w", err)
	}
	defer file.Close()

	var result [][]float64

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		parts := strings.Fields(line)

		// Skip empty lines
		if len(parts) == 0 {
			continue
		}

		// Validate number of columns
		if len(parts) != 2 {
			return nil, fmt.Errorf("invalid format in line: %q - expected 2 numbers, got %d", line, len(parts))
		}

		// Convert to float64
		x, err := strconv.ParseFloat(parts[0], 64)
		if err != nil {
			return nil, fmt.Errorf("error parsing float in line %q: %w", line, err)
		}

		y, err := strconv.ParseFloat(parts[1], 64)
		if err != nil {
			return nil, fmt.Errorf("error parsing float in line %q: %w", line, err)
		}

		result = append(result, []float64{x, y})
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading file: %w", err)
	}

	return result, nil
}

func getFilename(filePath string) string {
	// Get the base name (removes directory components)
	base := filepath.Base(filePath)

	// Remove the extension (everything after last dot)
	ext := filepath.Ext(base)

	// Trim the extension from base name
	nameWithoutExt := strings.TrimSuffix(base, ext)

	return nameWithoutExt
}

type CollisionEvent struct {
	x          int
	energyLoss float64
	collType   lxgata.CollisionType
	origin     int
}

type CSV [][]string

func (data CSV) Less(i, j int) bool {
	return natsort.Compare(data[i][0], data[j][0])
}

func (data CSV) Len() int {
	return len(data)
}
func (data CSV) Swap(i, j int) {
	data[i], data[j] = data[j], data[i]
}

func writeAsCSV(data CSV, path, subpath, filename string, columns []string) {
	var scalarW *csv.Writer
	println(filename, path)
	clearName := getFilename(filename)
	scalarParams, err := openFile(true, path, subpath, clearName)
	if err != nil {
		println("unable to save dc and secondary emission coefficient: ", err.Error())
		os.Exit(1)
	} else {
		scalarW = csv.NewWriter(scalarParams)
		scalarW.WriteAll([][]string{
			columns,
		})
		scalarW.Flush()
	}
	sort.Sort(data)
	scalarW.WriteAll(data)
	scalarW.Flush()
}
