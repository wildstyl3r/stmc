package main

import (
	"encoding/csv"
	"flag"
	"fmt"
	"log"
	"math"
	"os"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/BurntSushi/toml"
	"github.com/wildstyl3r/lxgata"
)

type output struct {
	saveFlag    *bool
	fileSuffix  string
	columnNames []string
	indexStep   *float64
	data        *[]float64
	scalers     []func(float64) float64
	notAverage  bool
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
	Runs                 int
	NElectrons           int
	MakeDir              bool
}

type Config struct {
	OutputDir string
	Models    map[string]ModelParameters
}

func main() {
	var TownsendAlpha []float64
	var TownsendAlphaD []float64
	var TownsendAlphaF []float64
	var driftVelocity []float64
	var electronDensity []float64
	var meanEnergy []float64
	var ionizationSourceTerm []float64
	var ionizationSourceTermD []float64
	var ionizations []float64
	var flux []float64
	var xStep float64
	var numCells int

	outputs := map[string]output{
		"Townsend alpha": {
			saveFlag:    flag.Bool("ta", false, "save Townsend alpha coefficient"),
			fileSuffix:  "Townsend_alpha",
			columnNames: []string{"px (Torr cm)", "a/p (cm^-1 Torr^-1)"},
			indexStep:   &xStep,
			data:        &TownsendAlpha,
			scalers:     []func(float64) float64{m2cm, cm2m},
		},
		"Townsend alpha per density": {
			saveFlag:    flag.Bool("tad", false, "save Townsend alpha coefficient divided by electron density"),
			fileSuffix:  "Townsend_alpha_D",
			columnNames: []string{"px (Torr cm)", "a*n/p (cm^2 Torr^-1)"},
			indexStep:   &xStep,
			data:        &TownsendAlphaD,
			scalers:     []func(float64) float64{m2cm, func(f float64) float64 { return f * 1.e4 }},
		},
		"Drift velocity": {
			saveFlag:    flag.Bool("vx", false, "save drift velocity"),
			fileSuffix:  "drift_vel",
			columnNames: []string{"px (Torr cm)", "v_x (cm s^-1)"},
			indexStep:   &xStep,
			data:        &driftVelocity,
			scalers:     []func(float64) float64{m2cm, m2cm},
		},
		"Density": {
			saveFlag:    flag.Bool("n", false, "save density"),
			fileSuffix:  "density",
			columnNames: []string{"px (Torr cm)", "n (сm^-3)"},
			indexStep:   &xStep,
			data:        &electronDensity,
			scalers:     []func(float64) float64{m2cm, func(f float64) float64 { return f / 1.e6 }},
		},
		"Mean energy": {
			saveFlag:    flag.Bool("emean", false, "save mean energy"),
			fileSuffix:  "mean_energy",
			columnNames: []string{"px (Torr cm)", "e (eV)"},
			indexStep:   &xStep,
			data:        &meanEnergy,
			scalers:     []func(float64) float64{m2cm},
		},
		"Ionization source term": {
			saveFlag:    flag.Bool("stf", false, "save ionization source term"),
			fileSuffix:  "ist",
			columnNames: []string{"px (Torr cm)", "S_n/p (cm^-1 Torr^-1)"},
			indexStep:   &xStep,
			data:        &ionizationSourceTerm,
			scalers:     []func(float64) float64{m2cm, cm2m},
		},
		"Ionization source term per density": {
			saveFlag:    flag.Bool("std", false, "save ionization source term divided by electron density"),
			fileSuffix:  "istd",
			columnNames: []string{"px (Torr cm)", "S_n*n/p (cm^2 Torr^-1)"},
			indexStep:   &xStep,
			data:        &ionizationSourceTermD,
			scalers:     []func(float64) float64{m2cm, func(f float64) float64 { return f * 1.e4 }},
		},
		"Ionization counter": {
			saveFlag:    flag.Bool("ic", false, "save ionization counter"),
			fileSuffix:  "ic",
			columnNames: []string{"px (Torr cm)", "N_i(cm^-1 Torr^-1)"},
			indexStep:   &xStep,
			data:        &ionizations,
			scalers:     []func(float64) float64{m2cm, cm2m},
		},
		"Flux": {
			saveFlag:    flag.Bool("f", false, "save flux"),
			fileSuffix:  "flux",
			columnNames: []string{"px (Torr cm)", "Phi(cm^-1 Torr^-1)"},
			indexStep:   &xStep,
			data:        &flux,
			scalers:     []func(float64) float64{m2cm, cm2m},
		},
		"Townsend Alpha from flux": {
			saveFlag:    flag.Bool("taf", false, "save Townsend alpha calculated as 1/F * dF/dx "),
			fileSuffix:  "Townsend_Alpha_F",
			columnNames: []string{"px (Torr cm)", "a/p (cm^-1 Torr^-1)"},
			indexStep:   &xStep,
			data:        &TownsendAlphaF,
			scalers:     []func(float64) float64{m2cm, cm2m},
			notAverage:  true,
		},
	}
	//var densityNorm = flag.Bool("dn", false, "divide TA by distribution integral")
	var configFileNamePointer = flag.String("input", "He_Tran_norm", "model configuration in toml format")
	flag.Parse()

	startTime := time.Now()
	fmt.Printf("Current time: %s\n", startTime.UTC().Format(time.UnixDate))

	configFileName := strings.TrimSuffix(*configFileNamePointer, ".toml")

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

	outputPath := ""

	if config.OutputDir != "" && config.OutputDir != "." {
		os.MkdirAll(config.OutputDir, 0750)
		outputPath += config.OutputDir + "/"
	}

	for modelName, parameters := range config.Models {
		fmt.Println("\n" + modelName)
		if !meta.IsDefined("Models", modelName, "CrossSections") ||
			!meta.IsDefined("Models", modelName, "GapLength") ||
			!meta.IsDefined("Models", modelName, "CathodeFallLength") ||
			!meta.IsDefined("Models", modelName, "CathodeFallPotential") {
			fmt.Println("Model " + modelName + " lacks key parameters (CrossSections, GapLength, CathodeFallLength or CathodeFallPotential)")
			continue
		}
		if !meta.IsDefined("Models", modelName, "CathodeCurrent") {
			parameters.CathodeCurrent = 2.84e-2
		}
		if !meta.IsDefined("Models", modelName, "ConstEField") {
			parameters.ConstEField = -100.
		}
		if !meta.IsDefined("Models", modelName, "Temperature") {
			parameters.Temperature = 300.
		}
		if !meta.IsDefined("Models", modelName, "DeltaE") {
			parameters.DeltaE = 0.3
		}
		if !meta.IsDefined("Models", modelName, "DeltaMu") {
			parameters.DeltaMu = 1. / 45.
		}
		if !meta.IsDefined("Models", modelName, "Runs") {
			parameters.Runs = 5
		}
		if !meta.IsDefined("Models", modelName, "NElectrons") {
			parameters.NElectrons = 100
		}
		if !meta.IsDefined("Models", modelName, "MakeDir") {
			parameters.MakeDir = false
		}

		crossSections, err := lxgata.LoadCrossSections(parameters.CrossSections)
		if err != nil {
			panic(fmt.Errorf("invalid cross section file: %w", err))
		}

		var chanWg sync.WaitGroup
		dataflow := make(chan *Model)
		for run := 0; run < parameters.Runs; run++ {
			chanWg.Add(1)
			//worker
			go func() {
				defer chanWg.Done()
				setup := Model{
					crossSections:     crossSections,
					cathodeFallLength: cm2m(parameters.CathodeFallLength),
					Vc:                -math.Abs(parameters.CathodeFallPotential),
					gapLength:         cm2m(parameters.GapLength),
					EConst:            parameters.ConstEField,
					density:           p / (k * parameters.Temperature),
					nElectrons:        int(parameters.NElectrons),
					eStep:             parameters.DeltaE,
					muStep:            parameters.DeltaMu,
					cathodeFlux:       parameters.CathodeCurrent / electronCharge,
				}
				setup.init()
				setup.run()
				dataflow <- &setup
			}()
		}

		// chan killer
		go func() {
			chanWg.Wait()
			close(dataflow)
		}()

		counter := 0
		fmt.Printf("\rDone:[0/%d]", parameters.Runs)
		first := true
		for results := range dataflow {
			fmt.Printf("\rDone:[%d/%d]", counter+1, parameters.Runs)
			if first {
				first = false
				xStep = results.xStep
				numCells = results.numCells

				TownsendAlpha = make([]float64, numCells)
				TownsendAlphaD = make([]float64, numCells)
				TownsendAlphaF = make([]float64, numCells)
				driftVelocity = make([]float64, numCells)
				meanEnergy = make([]float64, numCells)
				electronDensity = make([]float64, numCells)
				ionizationSourceTerm = make([]float64, numCells)
				ionizationSourceTermD = make([]float64, numCells)
				ionizations = make([]float64, numCells)
				flux = make([]float64, numCells)
			}

			for xIndex := 0; xIndex < numCells; xIndex++ {
				TownsendAlpha[xIndex] += results.TownsendAlpha[xIndex]
				ionizationSourceTerm[xIndex] += results.sourceTerm[xIndex]
				meanEnergy[xIndex] += results.electronEnergy[xIndex]
				driftVelocity[xIndex] += results.driftVelocity[xIndex]
				electronDensity[xIndex] += results.electronDensity[xIndex]
				ionizations[xIndex] += float64(results.ionizationAtCell[xIndex]) / float64(results.nElectrons)
				flux[xIndex] += float64(results.flux[xIndex]) / float64(results.nElectrons)
				TownsendAlphaD[xIndex] += results.TownsendAlpha[xIndex] / results.electronDensity[xIndex]
				ionizationSourceTermD[xIndex] += results.sourceTerm[xIndex] * (parameters.CathodeCurrent / electronCharge) / (results.electronDensity[xIndex] * results.driftVelocity[xIndex])
			}
			counter++
		}
		for xIndex := 0; xIndex < numCells; xIndex++ {
			for _, output := range outputs {
				if output.notAverage {
					continue
				}
				(*output.data)[xIndex] /= float64(counter)
			}
		}

		for xIndex := 1; xIndex+2 < numCells; xIndex++ {
			TownsendAlphaF[xIndex] = (flux[xIndex+1] + flux[xIndex+2] - flux[xIndex] - flux[xIndex-1]) / (2. * flux[xIndex] * xStep)
		}
		println()

		if parameters.MakeDir {
			err := os.MkdirAll(outputPath+modelName, 0750)
			if err != nil {
				fmt.Fprintln(os.Stderr, err)
				continue
			}
			modelName += "/"
		} else {
			modelName += "_"
		}

		for name, output := range outputs {
			if *output.saveFlag {
				file, err := os.Create(outputPath + modelName + output.fileSuffix + ".txt")
				if err != nil {
					println("unable to save "+name+": ", err)
				} else {
					rows := [][]string{output.columnNames}
					for index := range *output.data {
						xColumnValue := *output.indexStep * float64(index)
						yColumnValue := (*output.data)[index]
						if len(output.scalers) > 0 && output.scalers[0] != nil {
							xColumnValue = output.scalers[0](xColumnValue)
						}
						if len(output.scalers) > 1 && output.scalers[1] != nil {
							yColumnValue = output.scalers[1](yColumnValue)
						}

						rows = append(rows, []string{strconv.FormatFloat(xColumnValue, 'f', -1, 64), strconv.FormatFloat(yColumnValue, 'f', -1, 64)})
					}
					w := csv.NewWriter(file)
					w.WriteAll(rows)
					println(name + " saved")
					if err := w.Error(); err != nil {
						log.Fatalln("error writing csv:", err)
					}
				}
			}
		}
	}
	fmt.Printf("Elapsed time: %v\n", time.Since(startTime))
}
