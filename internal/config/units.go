package config

import (
	"fmt"
	"strconv"
	"strings"

	"github.com/wildstyl3r/stmc/internal/constants"
	"github.com/wildstyl3r/stmc/internal/utils"
)

var unitToSIeV = map[string]float64{
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
	"s":    1,
	"ms":   1e-3,
	"mks":  1e-6,
	"J":    1 / constants.ElementaryCharge,
	"eV":   1,
	"V":    1,
}

type UnitClass int

const (
	Length UnitClass = iota
	Current
	Pressure
	Energy
	Time
	Voltage
)

var string2class = map[string]UnitClass{
	"Pressure": Pressure,
	"Length":   Length,
	"Current":  Current,
	"Time":     Time,
	"Energy":   Energy,
	"Voltage":  Voltage,
}

var unitsInClass = map[UnitClass][]string{
	Length:   {"mm", "cm", "m"},
	Current:  {"mkA", "mA", "A"},
	Pressure: {"Torr", "mbar", "bar", "Pa"},
	Energy:   {"eV", "J"},
	Time:     {"mks", "ms", "s"},
	Voltage:  {"V"},
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
	"mks":  Time,
	"ms":   Time,
	"s":    Time,
	"eV":   Energy,
	"J":    Energy,
	"V":    Voltage,
}

type UnitElement struct {
	Class UnitClass
	Power int
}

func (ue UnitElement) String(units []string) string {
	for n := range units {
		if class, exist := classesOfUnits[units[n]]; exist && class == ue.Class && ue.Power != 0 {
			if ue.Power == 1 {
				return units[n]
			}
			return fmt.Sprintf("%s^%d", units[n], ue.Power)
		}
	}
	return ""
}

type Unit []UnitElement

func GetUnitElements(tag string) (u Unit) {
	classes := strings.Split(tag, ",")
	for c := range classes {
		description := strings.Split(classes[c], ":")
		if len(description) == 2 {
			power, err := strconv.Atoi(description[1])
			if err != nil {
				panic("Unit error: unable to parse power")
			}
			u = append(u, UnitElement{Class: string2class[description[0]], Power: power})
		}
	}
	return u
}

func (u Unit) String(units []string) (s string) {
	temp := []string{}
	for e := range u {
		temp = append(temp, u[e].String(units))
	}
	return "(" + strings.Join(temp, "*") + ")"
}

var defaultUnits = []string{"mkA", "cm", "Torr", "s", "eV", "V"}

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

func FieldToSIeV(v float64, classes Unit, units []string, direct bool) float64 {
	for i := range classes {
		uc := classes[i]
		unit := utils.AnyIntersection(unitsInClass[uc.Class], units)
		if unit == "" {
			continue
		}
		absPower := utils.IntAbs(uc.Power)
		if direct {
			if uc.Power > 0 {
				for range absPower {
					v *= unitToSIeV[unit]
				}
			} else {
				for range absPower {
					v /= unitToSIeV[unit]
				}
			}
		} else {
			if uc.Power > 0 {
				for range absPower {
					v /= unitToSIeV[unit]
				}
			} else {
				for range absPower {
					v *= unitToSIeV[unit]
				}
			}
		}
	}
	return v
}
