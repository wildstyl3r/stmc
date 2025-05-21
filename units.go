package main

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
	length UnitClass = iota
	current
	pressure
	energy
)

var unitsInClass = map[UnitClass][]string{
	length:   {"mm", "cm", "m"},
	current:  {"mkA", "mA", "A"},
	pressure: {"Torr", "mbar", "bar", "Pa"},
	energy:   {"eV", "J"},
}

var classesOfUnits = map[string]UnitClass{
	"Pa":   pressure,
	"bar":  pressure,
	"mbar": pressure,
	"Torr": pressure,
	"m":    length,
	"cm":   length,
	"mm":   length,
	"A":    current,
	"mA":   current,
	"mkA":  current,
}

type UnitElement = struct {
	class UnitClass
	power int
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
		unit := intersect(unitsInClass[uc.class], units)
		absPower := intAbs(uc.power)
		if direct {
			if uc.power > 0 {
				for range absPower {
					v *= unitToSI[*unit]
				}
			} else {
				for range absPower {
					v /= unitToSI[*unit]
				}
			}
		} else {
			if uc.power > 0 {
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
