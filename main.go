package main

import (
	"encoding/csv"
	"flag"
	"fmt"
	"log"
	"math"
	"os"
	"runtime/pprof"
	"strconv"
	"strings"
	"time"

	"github.com/BurntSushi/toml"
	"github.com/wildstyl3r/lxgata"
)

type output struct {
	saveFlag    *bool
	fileSuffix  string
	columnNames []string
	indexStep   *float64
	value       func(int) float64
	scalers     []func(float64) float64
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

//AINF
// AIV?

func main() {
	var cpuprofile = flag.String("cpuprofile", "", "write cpu profile to file")
	var TownsendAlphaF []float64
	var rateIntegral []float64
	var electronDensity []float64
	var meanEnergy []float64
	var ionizations []float64
	var flux []float64
	var xStep float64
	var Navg, Vxavg, Favg, Eavg float64
	var cathodeFlux float64
	var all = flag.Bool("all", false, "save every available metric")

	outputs := map[string]output{
		"Townsend alpha": {
			saveFlag:    flag.Bool("ta", false, "save Townsend alpha coefficient"),
			fileSuffix:  "Townsend_alpha",
			columnNames: []string{"px (Torr cm)", "a/p (cm^-1 Torr^-1)"},
			indexStep:   &xStep,
			value:       func(x int) float64 { return rateIntegral[x] * electronDensity[x] / flux[x] },
			scalers:     []func(float64) float64{m2cm, cm2m},
		},
		"Townsend alpha per density": {
			saveFlag:    flag.Bool("tad", false, "save Townsend alpha coefficient divided by electron density"),
			fileSuffix:  "Townsend_alpha_D",
			columnNames: []string{"px (Torr cm)", "a*n/p (cm^2 Torr^-1)"},
			indexStep:   &xStep,
			value:       func(x int) float64 { return rateIntegral[x] / flux[x] },
			scalers:     []func(float64) float64{m2cm, func(f float64) float64 { return f * 1.e4 }},
		},
		"Drift velocity": {
			saveFlag:    flag.Bool("vx", false, "save drift velocity"),
			fileSuffix:  "drift_vel",
			columnNames: []string{"px (Torr cm)", "v_x (cm s^-1)"},
			indexStep:   &xStep,
			value:       func(x int) float64 { return flux[x] / electronDensity[x] },
			scalers:     []func(float64) float64{m2cm, m2cm},
		},
		"Density": {
			saveFlag:    flag.Bool("n", false, "save density"),
			fileSuffix:  "density",
			columnNames: []string{"px (Torr cm)", "n (сm^-3)"},
			indexStep:   &xStep,
			value:       func(x int) float64 { return electronDensity[x] },
			scalers:     []func(float64) float64{m2cm, func(f float64) float64 { return f / 1.e6 }},
		},
		"Mean energy": {
			saveFlag:    flag.Bool("emean", false, "save mean energy"),
			fileSuffix:  "mean_energy",
			columnNames: []string{"px (Torr cm)", "e (eV)"},
			indexStep:   &xStep,
			value:       func(x int) float64 { return meanEnergy[x] },
			scalers:     []func(float64) float64{m2cm},
		},
		"Ionization source term": {
			saveFlag:    flag.Bool("stf", false, "save ionization source term"),
			fileSuffix:  "ist",
			columnNames: []string{"px (Torr cm)", "S_n/p (cm^-1 Torr^-1)"},
			indexStep:   &xStep,
			value:       func(x int) float64 { return rateIntegral[x] * electronDensity[x] / cathodeFlux },
			scalers:     []func(float64) float64{m2cm, cm2m},
		},
		"Ionization source term per density": {
			saveFlag:    flag.Bool("std", false, "save ionization source term divided by electron density"),
			fileSuffix:  "istd",
			columnNames: []string{"px (Torr cm)", "S_n*n/p (cm^2 Torr^-1)"},
			indexStep:   &xStep,
			value:       func(x int) float64 { return rateIntegral[x] / cathodeFlux },
			scalers:     []func(float64) float64{m2cm, func(f float64) float64 { return f * 1.e4 }},
		},
		"Ionization counter": {
			saveFlag:    flag.Bool("ic", false, "save ionization counter"),
			fileSuffix:  "ic",
			columnNames: []string{"px (Torr cm)", "N_i(cm^-1 Torr^-1)"},
			indexStep:   &xStep,
			value:       func(x int) float64 { return ionizations[x] },
			scalers:     []func(float64) float64{m2cm, cm2m},
		},
		"Electron flux": {
			saveFlag:    flag.Bool("f", false, "save calculated flux"),
			fileSuffix:  "flux",
			columnNames: []string{"px (Torr cm)", "Phi(cm^-1 Torr^-1)"},
			indexStep:   &xStep,
			value:       func(x int) float64 { return flux[x] },
			scalers:     []func(float64) float64{m2cm, cm2m},
		},
		// "Townsend Alpha from calculated flux": {
		// 	saveFlag:    flag.Bool("taf", false, "save Townsend alpha calculated as 1/(n(x)*v(x)) * d(n(x)v(x))/dx "),
		// 	fileSuffix:  "Townsend_Alpha_F",
		// 	columnNames: []string{"px (Torr cm)", "a/p (cm^-1 Torr^-1)"},
		// 	indexStep:   &xStep,
		// 	data:        &TownsendAlphaF,
		// 	scalers:     []func(float64) float64{m2cm, cm2m},
		// },
		"Ionization rate": {
			saveFlag:    flag.Bool("ir", false, "save ionization rate "),
			fileSuffix:  "ir",
			columnNames: []string{"px (Torr cm)", "nu (s^-1)"},
			indexStep:   &xStep,
			value:       func(x int) float64 { return rateIntegral[x] },
			scalers:     []func(float64) float64{m2cm},
		},
		"Ionization rate mul density": {
			saveFlag:    flag.Bool("ird", false, "save ionization rate mul density"),
			fileSuffix:  "irxd",
			columnNames: []string{"px (Torr cm)", "R (cm^-3 s^-1)"},
			indexStep:   &xStep,
			value:       func(x int) float64 { return rateIntegral[x] * electronDensity[x] },
			scalers:     []func(float64) float64{m2cm, func(f float64) float64 { return f / 1.e6 }},
		},
		"Ionization rate per flux": {
			saveFlag:    flag.Bool("irf", false, "save ionization rate per flux"),
			fileSuffix:  "irpf",
			columnNames: []string{"px (Torr cm)", "RpF (cm^2)"},
			indexStep:   &xStep,
			value:       func(x int) float64 { return rateIntegral[x] / flux[x] },
			scalers:     []func(float64) float64{m2cm, func(f float64) float64 { return f * 1.e4 }},
		},
		"Ionization rate per velocity": {
			saveFlag:    flag.Bool("irpv", false, "save ionization rate per velocity"),
			fileSuffix:  "irpv",
			columnNames: []string{"px (Torr cm)", "RpV (cm^-1)"},
			indexStep:   &xStep,
			value:       func(x int) float64 { return rateIntegral[x] * electronDensity[x] / flux[x] },
			scalers:     []func(float64) float64{m2cm, cm2m},
		},
		"Ionization rate per density": {
			saveFlag:    flag.Bool("irpd", false, "save ionization rate per density"),
			fileSuffix:  "irpd",
			columnNames: []string{"px (Torr cm)", "RpV (cm^3 s^-1)"},
			indexStep:   &xStep,
			value:       func(x int) float64 { return rateIntegral[x] / electronDensity[x] },
			scalers:     []func(float64) float64{m2cm, func(f float64) float64 { return f * 1.e6 }},
		},
		"I*Navg": {
			saveFlag:    flag.Bool("ain", false, "i*navg"),
			fileSuffix:  "ain",
			columnNames: []string{"px (Torr cm)", "??"},
			indexStep:   &xStep,
			value:       func(x int) float64 { return rateIntegral[x] * Navg },
			scalers:     []func(float64) float64{m2cm},
		},
		"I/Vavg": {
			saveFlag:    flag.Bool("aiv", false, "i/vavg"),
			fileSuffix:  "aiv",
			columnNames: []string{"px (Torr cm)", "??"},
			indexStep:   &xStep,
			value:       func(x int) float64 { return rateIntegral[x] / Vxavg },
			scalers:     []func(float64) float64{m2cm},
		},
		"I/Favg": {
			saveFlag:    flag.Bool("aif", false, "i/favg"),
			fileSuffix:  "aif",
			columnNames: []string{"px (Torr cm)", "??"},
			indexStep:   &xStep,
			value:       func(x int) float64 { return rateIntegral[x] / Favg },
			scalers:     []func(float64) float64{m2cm},
		},
		"I*Navg/Favg": {
			saveFlag:    flag.Bool("ainf", false, "i*navg"),
			fileSuffix:  "ainf",
			columnNames: []string{"px (Torr cm)", "??"},
			indexStep:   &xStep,
			value:       func(x int) float64 { return rateIntegral[x] * Navg / Favg },
			scalers:     []func(float64) float64{m2cm},
		},
		"I/VavgFx": {
			saveFlag:    flag.Bool("aivx", false, "i/vavg*f(x)"),
			fileSuffix:  "aivx",
			columnNames: []string{"px (Torr cm)", "??"},
			indexStep:   &xStep,
			value:       func(x int) float64 { return rateIntegral[x] / Vxavg * flux[x] },
			scalers:     []func(float64) float64{m2cm},
		},
		"I*Navg/FavgFx": {
			saveFlag:    flag.Bool("ainfx", false, "i*navg/favg*f(x)"),
			fileSuffix:  "ainfx",
			columnNames: []string{"px (Torr cm)", "??"},
			indexStep:   &xStep,
			value:       func(x int) float64 { return rateIntegral[x] * Navg / Favg * flux[x] },
			scalers:     []func(float64) float64{m2cm},
		},
	}
	var configFileNamePointer = flag.String("input", "He_Tran_norm", "model configuration in toml format")
	flag.Parse()
	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

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
		noParams := false
		if !meta.IsDefined("Models", modelName, "CrossSections") {
			if meta.IsDefined("CrossSections") {
				parameters.CrossSections = config.CrossSections
			} else {
				noParams = true
			}
		}
		if !meta.IsDefined("Models", modelName, "GapLength") {
			if meta.IsDefined("GapLength") {
				parameters.GapLength = config.GapLength
			} else {
				noParams = true
			}
		}
		if !meta.IsDefined("Models", modelName, "CathodeFallLength") {
			if meta.IsDefined("CathodeFallLength") {
				parameters.CathodeFallLength = config.CathodeFallLength
			} else {
				noParams = true
			}
		}
		if !meta.IsDefined("Models", modelName, "CathodeFallPotential") {
			if meta.IsDefined("CathodeFallPotential") {
				parameters.CathodeFallPotential = config.CathodeFallPotential
			} else {
				noParams = true
			}
		}

		if noParams {
			fmt.Println("Model " + modelName + " lacks key parameters (CrossSections, GapLength, CathodeFallLength or CathodeFallPotential)")
			continue
		}

		if !meta.IsDefined("Models", modelName, "CathodeCurrent") {
			if meta.IsDefined("CathodeCurrent") {
				parameters.CathodeCurrent = config.CathodeCurrent
			} else {
				parameters.CathodeCurrent = 2.84e-2
			}
		}
		if !meta.IsDefined("Models", modelName, "ConstEField") {
			if meta.IsDefined("ConstEField") {
				parameters.ConstEField = config.ConstEField
			} else {
				parameters.ConstEField = -100.
			}
		}
		if !meta.IsDefined("Models", modelName, "Temperature") {
			if meta.IsDefined("Temperature") {
				parameters.Temperature = config.Temperature
			} else {
				parameters.Temperature = 300.
			}
		}
		if !meta.IsDefined("Models", modelName, "DeltaE") {
			if meta.IsDefined("DeltaE") {
				parameters.DeltaE = config.DeltaE
			} else {
				parameters.DeltaE = 0.3
			}
		}
		if !meta.IsDefined("Models", modelName, "DeltaMu") {
			if meta.IsDefined("DeltaMu") {
				parameters.DeltaMu = config.DeltaMu
			} else {
				parameters.DeltaMu = 1. / 45.
			}
		}
		if !meta.IsDefined("Models", modelName, "NElectrons") {
			if meta.IsDefined("NElectrons") {
				parameters.NElectrons = config.NElectrons
			} else {
				parameters.NElectrons = 100
			}
		}
		if !meta.IsDefined("Models", modelName, "MakeDir") {
			parameters.MakeDir = false
			if meta.IsDefined("MakeDir") {
				parameters.MakeDir = config.MakeDir
			} else {
				parameters.MakeDir = false
			}
		}

		crossSections, err := lxgata.LoadCrossSections(parameters.CrossSections)
		if err != nil {
			panic(fmt.Errorf("invalid cross section file: %w", err))
		}

		cathodeFlux = parameters.CathodeCurrent / electronCharge

		model := Model{
			crossSections:     crossSections,
			cathodeFallLength: cm2m(parameters.CathodeFallLength),
			Vc:                -math.Abs(parameters.CathodeFallPotential),
			gapLength:         cm2m(parameters.GapLength),
			EConst:            parameters.ConstEField,
			density:           p / (k * parameters.Temperature),
			nElectrons:        int(parameters.NElectrons),
			eStep:             parameters.DeltaE,
			muStep:            parameters.DeltaMu,
			cathodeFlux:       cathodeFlux,
		}
		model.init()
		model.run()

		electronDensity = make([]float64, model.numCells)
		rateIntegral = make([]float64, model.numCells)
		meanEnergy = make([]float64, model.numCells)
		flux = make([]float64, model.numCells)

		for xIndex := 0; xIndex < model.numCells; xIndex++ {
			for eIndex := 1; eIndex < model.numCellsE; eIndex++ {
				currentEnergy := model.eStep * float64(eIndex)
				fXE := 0.
				for muIndex := 0; muIndex < model.numCellsMu; muIndex++ {
					currentMu := model.muStep*float64(muIndex) - 1.

					f := 0.
					if currentMu > 0.0001 {
						f = model.psiFIncrement * float64(model.distribution[xIndex][eIndex][muIndex]) / (model.lookUpVelocity[eIndex] * currentMu)
					}

					electronDensity[xIndex] += f

					fXE += f // * s.muStep

					if eIndex > 0 {
						meanEnergy[xIndex] += f * currentEnergy

						flux[xIndex] += f * model.lookUpVelocity[eIndex] * currentMu
						Favg += flux[xIndex]

					}
				}
				cs := model.crossSections.TotalCrossSectionOfKindAt(lxgata.IONIZATION, currentEnergy)
				vel := model.lookUpVelocity[eIndex] //math.Sqrt(eV2J(currentEnergy)) * energyRoot2Velocity
				rateIntegral[xIndex] += cs * vel * fXE
			}
			meanEnergy[xIndex] /= electronDensity[xIndex]
			Vxavg += flux[xIndex] / electronDensity[xIndex]
			Navg += electronDensity[xIndex]
			Eavg += meanEnergy[xIndex]

		}
		Vxavg /= float64(model.numCells)
		Navg /= float64(model.numCells)
		Favg /= float64(model.numCells)
		Eavg /= float64(model.numCells)
		// constN := model.cathodeFlux / Vxavg
		// initV := math.Sqrt(2. * eV2J(4.5) / me)

		xStep = model.xStep

		TownsendAlphaF = make([]float64, model.numCells)
		ionizations = make([]float64, model.numCells)

		initElectronDensity := model.cathodeFlux / math.Sqrt(2.*eV2J(4.5)/me)
		fmt.Printf("initN: %v\n", initElectronDensity)
		for xIndex := 0; xIndex < model.numCells; xIndex++ {
			ionizations[xIndex] = float64(model.ionizationAtCell[xIndex]) / float64(model.nElectrons) * initElectronDensity
		}

		for xIndex := 1; xIndex+2 < model.numCells; xIndex++ {
			TownsendAlphaF[xIndex] = (flux[xIndex+1] + flux[xIndex+2] - flux[xIndex] - flux[xIndex-1]) / (2. * flux[xIndex] * xStep)
		}

		for name, output := range outputs {
			if *output.saveFlag || *all {
				var file *os.File
				var err error
				if parameters.MakeDir && output.fileSuffix != "" && output.fileSuffix != "." {
					os.MkdirAll(outputPath+output.fileSuffix, 0750)
					file, err = os.Create(outputPath + output.fileSuffix + "/" + modelName + ".txt")
				} else {
					file, err = os.Create(outputPath + modelName + "_" + output.fileSuffix + ".txt")
				}
				if err != nil {
					println("unable to save "+name+": ", err)
				} else {
					rows := [][]string{output.columnNames}
					for index := 0; index < model.numCells; index++ {
						xColumnValue := *output.indexStep * float64(index)
						yColumnValue := output.value(index)
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
