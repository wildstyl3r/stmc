package config

import "github.com/wildstyl3r/stmc/internal/utils"

var unitToSI = map[string]float64{
	"Pa":   1,              // [Pa]
	"bar":  1e5,            // [Pa]
	"mbar": 1e2,            // [Pa]
	"Torr": 101325. / 760., // [Pa]
	"m":    1,              // [m]
	"cm":   1e-2,           // [m]
	"mm":   1e-3,           // [m]
	"A":    1,              // [A]
	"mA":   1e-3,           // [A]
	"mkA":  1e-6,           // [A]
}

type UnitClass int

const (
	Length UnitClass = iota
	Current
	Pressure
	Energy
)

var unitsInClass = map[UnitClass][]string{
	Length:   {"mm", "cm", "m"},
	Current:  {"mkA", "mA", "A"},
	Pressure: {"Torr", "mbar", "bar", "Pa"},
	Energy:   {"eV", "J"},
}

var classesOfUnits = map[string]UnitClass{
	"Pa":   Pressure,
	"bar":  Pressure,
	"mbar": Pressure,
	"Torr": Pressure,
	"m":    Length,
	"cm":   Length,
	"mm":   Length,
	"A":    Current,
	"mA":   Current,
	"mkA":  Current,
}

type UnitElement = struct {
	Class UnitClass
	Power int
}

func checkUnits(units []string) (extended, conflicts []string) {
	classes := map[UnitClass]struct{}{}
	for _, unit := range units {
		if _, some := classes[classesOfUnits[unit]]; some {
			conflicts = append(conflicts, unit)
		} else {
			classes[classesOfUnits[unit]] = struct{}{}
		}
	}
	extended = units
	for _, unit := range defaultUnits {
		if _, some := classes[classesOfUnits[unit]]; !some {
			extended = append(extended, unit)
		}
	}
	return
}

func SI(v float64, classes []UnitElement, units []string, direct bool) float64 {
	for i := range classes {
		uc := classes[i]
		unit := utils.Intersect(unitsInClass[uc.Class], units)
		absPower := utils.IntAbs(uc.Power)
		if direct {
			if uc.Power > 0 {
				for range absPower {
					v *= unitToSI[*unit]
				}
			} else {
				for range absPower {
					v /= unitToSI[*unit]
				}
			}
		} else {
			if uc.Power > 0 {
				for range absPower {
					v /= unitToSI[*unit]
				}
			} else {
				for range absPower {
					v *= unitToSI[*unit]
				}
			}
		}
	}
	return v
}
