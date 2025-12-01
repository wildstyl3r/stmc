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

func (ue UnitElement) String(units UnitConfig) string {
	unit := units.GetForClass(ue.Class)
	if unit != "" {
		if ue.Power == 1 {
			return unit
		}
		return fmt.Sprintf("%s^%d", unit, ue.Power)
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

func (u Unit) String(units UnitConfig) (s string) {
	temp := []string{}
	for e := range u {
		temp = append(temp, u[e].String(units))
	}
	return "(" + strings.Join(temp, "*") + ")"
}

var defaultUnits = UnitConfig{
	L: "cm",
	I: "mkA",
	P: "Torr",
	T: "s",
	E: "eV",
} // []string{"mkA", "cm", "Torr", "s", "eV", "V"}

func mergeEmpty(target, source UnitConfig) (result UnitConfig) {
	if target.L == "" {
		target.L = source.L
	}
	if target.I == "" {
		target.I = source.I
	}
	if target.P == "" {
		target.P = source.P
	}
	if target.E == "" {
		target.E = source.E
	}
	if target.T == "" {
		target.T = source.E
	}
	return target
}

func checkUnits(units UnitConfig) (extended UnitConfig, unknowns []string) {
	units = mergeEmpty(units, defaultUnits)

	if _, known := classesOfUnits[units.L]; !known {
		unknowns = append(unknowns, units.L)
	}
	if _, known := classesOfUnits[units.I]; !known {
		unknowns = append(unknowns, units.I)
	}
	if _, known := classesOfUnits[units.I]; !known {
		unknowns = append(unknowns, units.I)
	}
	if _, known := classesOfUnits[units.T]; !known {
		unknowns = append(unknowns, units.T)
	}
	if _, known := classesOfUnits[units.P]; !known {
		unknowns = append(unknowns, units.P)
	}
	return units, unknowns
}

func FieldToSIeV(v float64, classes Unit, units UnitConfig, direct bool) float64 {
	for i := range classes {
		uc := classes[i]
		unit := units.GetForClass(uc.Class)
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
