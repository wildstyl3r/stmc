package main

import (
	"fmt"
	"os"

	"github.com/BurntSushi/toml"
)

type Config struct {
	OutputDir string
	Models    map[string]ModelParameters

	// to reset global defaults
	CrossSections        string
	GapLength            float64 // [m]
	CathodeFallLength    float64 // [m]
	CathodeFallPotential float64 // [V]
	CathodeCurrent       float64 // [A m^-2] = [C s^-1 m^-2]
	ConstEField          float64 // [V / m]
	Temperature          float64 // [K]
	DeltaE               float64 // [eV]
	DeltaMu              float64
	NElectrons           int
	MakeDir              bool
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
	CrossSections        string
	GapLength            float64 // [m]
	CathodeFallLength    float64 // [m]
	CathodeFallPotential float64 // [V]
	CathodeCurrent       float64 // [A m^-2] = [C s^-1 m^-2]
	ConstEField          float64 // [V / m]
	Temperature          float64 // [K]
	DeltaE               float64 // [eV]
	DeltaMu              float64
	NElectrons           int
	MakeDir              bool
}

func (mp *ModelParameters) checkDefaults(modelName string, config *Config, meta *toml.MetaData) bool {
	noParams := false
	if !meta.IsDefined("Models", modelName, "CrossSections") {
		if meta.IsDefined("CrossSections") {
			mp.CrossSections = config.CrossSections
		} else {
			noParams = true
		}
	}
	if !meta.IsDefined("Models", modelName, "GapLength") {
		if meta.IsDefined("GapLength") {
			mp.GapLength = config.GapLength
		} else {
			noParams = true
		}
	}
	if !meta.IsDefined("Models", modelName, "CathodeFallLength") {
		if meta.IsDefined("CathodeFallLength") {
			mp.CathodeFallLength = config.CathodeFallLength
		} else {
			noParams = true
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

	if !meta.IsDefined("Models", modelName, "CathodeCurrent") {
		if meta.IsDefined("CathodeCurrent") {
			mp.CathodeCurrent = config.CathodeCurrent
		} else {
			mp.CathodeCurrent = 2.84e-2
		}
	}
	if !meta.IsDefined("Models", modelName, "ConstEField") {
		if meta.IsDefined("ConstEField") {
			mp.ConstEField = config.ConstEField
		} else {
			mp.ConstEField = -100.
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
			mp.DeltaE = 0.3
		}
	}
	if !meta.IsDefined("Models", modelName, "DeltaMu") {
		if meta.IsDefined("DeltaMu") {
			mp.DeltaMu = config.DeltaMu
		} else {
			mp.DeltaMu = 1. / 45.
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
	return true
}
