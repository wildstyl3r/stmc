package main

import (
	"fmt"
	"os"

	"github.com/BurntSushi/toml"
)

type Config struct {
	OutputDir     string
	Models        map[string]ModelParameters
	CrossSections string

	// to reset global defaults
	ModelParameters
	// CathodeDistance            float64 // [сm]
}

func loadConfig(configFileName string) (Config, toml.MetaData) {
	var config Config
	meta, err := toml.DecodeFile(configFileName+".toml", &config)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}

	if len(config.Models) == 0 {
		fmt.Println("No models provided")
		os.Exit(0)
	}
	return config, meta
}

type ModelParameters struct {
	GapLength                        float64 // [сm]
	CathodeFallLength                float64 // [сm]
	CathodeFallPotential             float64 // [V]
	CathodeCurrentDensity            float64 // [mqA cm^-2] // ==>[A m^-2] == [C s^-1 m^-2]
	CathodeCurrentPerPressureSquared float64 // [mA/mbar^2]
	ConstEField                      float64 // [V / m]
	Temperature                      float64 // [K]
	Pressure                         float64 // [Pa]
	DAmbipolar                       float64
	AmbipolarCharacterScale          float64
	CathodeArea                      float64 // [cm^2]

	DeltaE                     float64 // [eV]
	DeltaMu                    float64
	NElectrons                 int
	MakeDir                    bool
	ParallelPlaneHollowCathode bool

	PressureInmBars bool
	PressureInTorrs bool
	SI              bool

	CalculateCathodeFallLength bool
	CathodeFallLengthPrecision float64
	LR                         float64
}

func (mp *ModelParameters) toSI() {
	if !mp.SI {
		mp.GapLength *= 0.01
		mp.CathodeFallLength *= 0.01

		mp.CathodeCurrentDensity *= (1e-6 * 1e2 * 1e2)
		mp.ConstEField *= 1.e2

		if mp.PressureInmBars {
			mp.Pressure *= 1e-3 * 1e5
		} else if mp.PressureInTorrs {
			mp.Pressure *= Torr
		}
	}
}

func (mp *ModelParameters) checkMetrizeSetDefaults(modelName string, config *Config, meta *toml.MetaData) bool {
	noParams := false
	if !meta.IsDefined("CrossSections") {
		noParams = true
	}
	if !meta.IsDefined("Models", modelName, "GapLength") {
		if meta.IsDefined("GapLength") {
			mp.GapLength = config.GapLength
		} else {
			noParams = true
		}
	}
	if !meta.IsDefined("Models", modelName, "CalculateCathodeFallLength") {
		if meta.IsDefined("CalculateCathodeFallLength") {
			mp.CalculateCathodeFallLength = config.CalculateCathodeFallLength
		}
	}

	if !meta.IsDefined("Models", modelName, "LR") {
		if meta.IsDefined("LR") {
			mp.LR = config.LR
		} else {
			mp.LR = 0.01
		}
	}

	if !meta.IsDefined("Models", modelName, "CathodeFallLength") {
		if meta.IsDefined("CathodeFallLength") {
			mp.CathodeFallLength = config.CathodeFallLength
		} else if !mp.CalculateCathodeFallLength {
			noParams = true
			// if meta.IsDefined("Gamma") && meta.IsDefined("CathodeCurrentDensity") {
			// 	mp.CathodeFallLength = math.Sqrt(2.*(1.+config.Gamma)*driftVelHelium(ENRatio float64)*mp.CathodeFallPotential*e0 / config.CathodeCurrentDensity)
			// } else {
			// 	noParams = true
			// }
		} else {
			if !meta.IsDefined("Models", modelName, "CathodeFallLengthPrecision") {
				if meta.IsDefined("CathodeFallLengthPrecision") {
					mp.CathodeFallLengthPrecision = config.CathodeFallLengthPrecision
				} else {
					mp.CathodeFallLengthPrecision = 5e-6
				}
			}
		}
	}
	if !meta.IsDefined("Models", modelName, "CathodeFallPotential") {
		if meta.IsDefined("CathodeFallPotential") {
			mp.CathodeFallPotential = config.CathodeFallPotential
		} else {
			noParams = true
		}
	}

	if noParams {
		fmt.Println("Model " + modelName + " lacks key parameters (CrossSections, GapLength, CathodeFallLength or CathodeFallPotential)")
		return false
	}

	if !meta.IsDefined("Models", modelName, "CathodeCurrentDensity") {
		if meta.IsDefined("CathodeCurrentDensity") {
			mp.CathodeCurrentDensity = config.CathodeCurrentDensity
		} else {
			mp.CathodeCurrentDensity = 2.84 //[mqA / cm^2]
		}
	}

	if !meta.IsDefined("Models", modelName, "Pressure") {
		if meta.IsDefined("Pressure") {
			mp.Pressure = config.Pressure
		} else {
			mp.Pressure = 101325. / 760. // this is 133.322 Pa = 1 Torr
		}
	}

	if meta.IsDefined("Models", modelName, "CathodeCurrentPerPressureSquared") {
		if meta.IsDefined("Models", modelName, "CathodeArea") {
			//mA*mbar^-2                            mbar         mbar            cm^-2
			mp.CathodeCurrentDensity = 1000. * mp.CathodeCurrentPerPressureSquared * mp.Pressure * mp.Pressure / mp.CathodeArea
			// * 1000 to translate from [mA * cm^-2] to [mqA * cm^-2]
		} else if meta.IsDefined("CathodeArea") {
			mp.CathodeCurrentDensity = 1000. * mp.CathodeCurrentPerPressureSquared * mp.Pressure * mp.Pressure / config.CathodeArea
		}

	}

	if !meta.IsDefined("Models", modelName, "ConstEField") {
		if meta.IsDefined("ConstEField") {
			mp.ConstEField = config.ConstEField
		} else {
			mp.ConstEField = -1.
		}
	}
	if !meta.IsDefined("Models", modelName, "Temperature") {
		if meta.IsDefined("Temperature") {
			mp.Temperature = config.Temperature
		} else {
			mp.Temperature = 300.
		}
	}
	if !meta.IsDefined("Models", modelName, "DeltaE") {
		if meta.IsDefined("DeltaE") {
			mp.DeltaE = config.DeltaE
		} else {
			mp.DeltaE = 0.1
		}
	}
	if !meta.IsDefined("Models", modelName, "DeltaMu") {
		if meta.IsDefined("DeltaMu") {
			mp.DeltaMu = config.DeltaMu
		} else {
			mp.DeltaMu = 1. / 90.
		}
	}
	if !meta.IsDefined("Models", modelName, "NElectrons") {
		if meta.IsDefined("NElectrons") {
			mp.NElectrons = config.NElectrons
		} else {
			mp.NElectrons = 100
		}
	}
	if !meta.IsDefined("Models", modelName, "MakeDir") {
		mp.MakeDir = true
		if meta.IsDefined("MakeDir") {
			mp.MakeDir = config.MakeDir
		}
	}
	if !meta.IsDefined("Models", modelName, "ParallelPlaneHollowCathode") {
		mp.ParallelPlaneHollowCathode = false
		if meta.IsDefined("ParallelPlaneHollowCathode") {
			mp.ParallelPlaneHollowCathode = config.ParallelPlaneHollowCathode
		}
	}

	if mp.ParallelPlaneHollowCathode {
		mp.GapLength /= 2
	}

	if !meta.IsDefined("Models", modelName, "PressureInmBars") {
		mp.PressureInmBars = false
		if meta.IsDefined("PressureInmBars") {
			mp.PressureInmBars = config.PressureInmBars
		}
	}
	println("inmbars", config.PressureInmBars, mp.PressureInmBars)
	if !meta.IsDefined("Models", modelName, "PressureInTorrs") {
		mp.PressureInTorrs = false
		if meta.IsDefined("PressureInTorrs") {
			mp.PressureInTorrs = config.PressureInTorrs
		}
	}
	mp.toSI()
	return true
}
