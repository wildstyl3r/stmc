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

	"github.com/spf13/viper"
	"github.com/wildstyl3r/lxgata"
)

func main() {
	// if len(os.Args) < 2 {
	// 	panic(fmt.Errorf("no input provided"))
	// }

	var saveTownsendAlpha = flag.Bool("ta", false, "save Townsend alpha coefficient")
	var saveDriftVelocity = flag.Bool("vx", false, "save drift velocity")
	var saveElectronDensity = flag.Bool("n", false, "save electron density")
	var saveEnergy = flag.Bool("emean", false, "save mean energy")
	var saveSourceTermFormula = flag.Bool("stf", false, "save source term calculated by paper formula")
	var configFileNamePointer = flag.String("input", "He_Tran_norm", "model configuration in toml format")
	flag.Parse()
	// var saveSourceTermFormula = flag.Bool("stf", false, "save source term calculated by paper formula")
	// var saveSourceTermNative = flag.Bool("stn", false, "save source term counted natively")

	configFileName := strings.TrimSuffix(*configFileNamePointer, ".toml")

	viper.SetConfigName(configFileName)
	viper.SetConfigType("toml") // REQUIRED if the config file does not have the extension in the name
	viper.AddConfigPath(".")    // optionally look for config in the working directory

	viper.SetDefault("temperature", 300.)
	viper.SetDefault("num_runs", 5)
	viper.SetDefault("deltaE", 0.3) // eV
	viper.SetDefault("deltaMu", 2./90.)
	viper.SetDefault("n_electrons", 100)
	viper.SetDefault("const_E_field", -100.)
	viper.SetDefault("cathode_current", 2.84e-2) // [A m^-2] = [C / s^-1 m^-2]

	err := viper.ReadInConfig() // Find and read the config file
	if err != nil {             // Handle errors reading the config file
		panic(fmt.Errorf("fatal error config file: %w", err))
	}

	//var angularDistributionAveraged [][][]float64
	var TownsendAlpha []float64
	var driftVelocity []float64
	var electronDensity []float64
	var meanEnergy []float64
	var normalizedIonizationSourceTerm []float64
	var xStep float64

	numRuns := viper.GetInt("num_runs")

	crossSections, err := lxgata.LoadCrossSections(viper.GetString("cross_sections"))
	if err != nil {
		panic(fmt.Errorf("invalid cross section file: %w", err))
	}

	startTime := time.Now()
	fmt.Printf("Current time: %s\n", startTime.UTC().Format(time.UnixDate))

	var chanWg sync.WaitGroup
	dataflow := make(chan *Model)
	for run := 0; run < numRuns; run++ {
		chanWg.Add(1)
		//worker
		go func() {
			defer chanWg.Done()
			setup := Model{
				crossSections:     crossSections,
				cathodeFallLength: cm2m(viper.GetFloat64("cathode_fall_length")),
				Vc:                -math.Abs(viper.GetFloat64("cathode_fall_potential")),
				gapLength:         cm2m(viper.GetFloat64("gap_length")),
				EConst:            viper.GetFloat64("const_E_field"),
				density:           p / (k * viper.GetFloat64("temperature")),
				nElectrons:        int(viper.GetFloat64("n_electrons")),
				eStep:             viper.GetFloat64("deltaE"),
				muStep:            viper.GetFloat64("deltaMu"),
				cathodeFlux:       viper.GetFloat64("cathode_current") / electronCharge,
			}
			setup.init()
			setup.run()
			dataflow <- &setup
		}()
	}

	// chanel killer
	go func() {
		chanWg.Wait()
		close(dataflow)
	}()

	// merger
	//go func() {
	counter := 0
	fmt.Printf("\rDone:[0/%d]", numRuns)
	for results := range dataflow {
		fmt.Printf("\rDone:[%d/%d]", counter+1, numRuns)
		if meanEnergy == nil {
			xStep = results.xStep
			//eStep = setup.eStep
			//muStep = setup.muStep

			TownsendAlpha = make([]float64, results.numCells)
			driftVelocity = make([]float64, results.numCells)
			meanEnergy = make([]float64, results.numCells)
			electronDensity = make([]float64, results.numCells)
			normalizedIonizationSourceTerm = make([]float64, results.numCells)
		}

		for xIndex := 0; xIndex < results.numCells; xIndex++ {
			TownsendAlpha[xIndex] += results.TownsendAlpha[xIndex]
			normalizedIonizationSourceTerm[xIndex] += results.TownsendAlpha[xIndex] * float64(results.flux[xIndex])
			meanEnergy[xIndex] += results.electronEnergy[xIndex]
			driftVelocity[xIndex] += results.driftVelocity[xIndex]
			electronDensity[xIndex] += results.electronDensity[xIndex]
		}
		counter++
	}
	for xIndex := 0; xIndex < len(meanEnergy); xIndex++ {
		TownsendAlpha[xIndex] /= float64(counter)
		normalizedIonizationSourceTerm[xIndex] /= float64(counter)
		driftVelocity[xIndex] /= float64(counter)
		meanEnergy[xIndex] /= float64(counter)
		electronDensity[xIndex] /= float64(counter)
	}
	println()
	//}()

	// for run := 0; run < numRuns; run++ {
	// 	fmt.Println("run: ", run)

	// 	setup := Model{
	// 		crossSections:     crossSections,
	// 		cathodeFallLength: cm2m(viper.GetFloat64("cathode_fall_length")),
	// 		Vc:                -math.Abs(viper.GetFloat64("cathode_fall_potential")),
	// 		gapLength:         cm2m(viper.GetFloat64("gap_length")),
	// 		EConst:            viper.GetFloat64("const_E_field"),
	// 		density:           p / (k * viper.GetFloat64("temperature")),
	// 		nElectrons:       int(viper.GetFloat64("n_electrons")),
	// 		eStep:             viper.GetFloat64("deltaE"),
	// 		muStep:            viper.GetFloat64("deltaMu"),
	// 	}
	// 	setup.init()

	// 	if meanEnergy == nil {
	// 		fmt.Println("parameters loaded, variables initialized")
	// 		xStep = setup.xStep
	// 		//eStep = setup.eStep
	// 		//muStep = setup.muStep

	// 		TownsendAlpha = make([]float64, setup.numCells)
	// 		normalizedIonizationSourceTerm = make([]float64, setup.numCells)
	// 		driftVelocity = make([]float64, setup.numCells)
	// 		meanEnergy = make([]float64, setup.numCells)
	// 		electronDensity = make([]float64, setup.numCells)
	// 	}
	// 	setup.run()

	// 	for xIndex := 0; xIndex < setup.numCells; xIndex++ {
	// 		TownsendAlpha[xIndex] += setup.TownsendAlpha[xIndex]
	// 		normalizedIonizationSourceTerm[xIndex] += setup.TownsendAlpha[xIndex] * float64(setup.flux[xIndex])
	// 		meanEnergy[xIndex] += setup.electronEnergy[xIndex]
	// 		driftVelocity[xIndex] += setup.driftVelocity[xIndex]
	// 		electronDensity[xIndex] += setup.electronDensity[xIndex]
	// 	}
	// }
	// for xIndex := 0; xIndex < len(meanEnergy); xIndex++ {
	// 	TownsendAlpha[xIndex] /= float64(numRuns)
	// 	normalizedIonizationSourceTerm[xIndex] /= float64(numRuns)
	// 	driftVelocity[xIndex] /= float64(numRuns)
	// 	meanEnergy[xIndex] /= float64(numRuns)
	// 	electronDensity[xIndex] /= float64(numRuns)
	// }

	if *saveTownsendAlpha {
		file, err := os.Create(configFileName + "_Townsend_alpha.txt")
		if err != nil {
			println("unable to save first Townsend coefficient: ", err)
		} else {
			output := [][]string{{"px (Torr cm)", "a/p (cm^-1 Torr^-1)"}}
			for xIndex := range TownsendAlpha {
				output = append(output, []string{strconv.FormatFloat(m2cm(xStep*float64(xIndex)), 'f', -1, 64), strconv.FormatFloat(TownsendAlpha[xIndex]/100., 'f', -1, 64)})
			}
			w := csv.NewWriter(file)
			w.WriteAll(output)
			println("TA saved")
			if err := w.Error(); err != nil {
				log.Fatalln("error writing csv:", err)
			}
		}
	}

	if *saveDriftVelocity {
		file, err := os.Create(configFileName + "_drift_vel.txt")
		if err != nil {
			println("unable to save drift velocity: ", err)
		} else {
			output := [][]string{{"px (Torr cm)", "v_x (cm s^-1)"}}
			for xIndex := range driftVelocity {
				output = append(output, []string{strconv.FormatFloat(m2cm(xStep*float64(xIndex)), 'f', -1, 64), strconv.FormatFloat(m2cm(driftVelocity[xIndex]), 'f', -1, 64)})
			}
			w := csv.NewWriter(file)
			w.WriteAll(output)
			println("drift velocity saved")
			if err := w.Error(); err != nil {
				log.Fatalln("error writing csv:", err)
			}
		}
	}

	if *saveEnergy {
		file, err := os.Create(configFileName + "_mean_energy.txt")
		if err != nil {
			println("unable to save drift velocity: ", err)
		} else {
			output := [][]string{{"px (Torr cm)", "e (eV)"}}
			for xIndex := range meanEnergy {
				output = append(output, []string{strconv.FormatFloat(m2cm(xStep*float64(xIndex)), 'f', -1, 64), strconv.FormatFloat(meanEnergy[xIndex], 'f', -1, 64)})
			}
			w := csv.NewWriter(file)
			w.WriteAll(output)
			println("Energy saved")
			if err := w.Error(); err != nil {
				log.Fatalln("error writing csv:", err)
			}
		}
	}

	if *saveElectronDensity {
		file, err := os.Create(configFileName + "_density.txt")
		if err != nil {
			println("unable to save density: ", err)
		} else {
			output := [][]string{{"px (Torr cm)", "n (сm^-3)"}}
			for xIndex := range meanEnergy {
				output = append(output, []string{strconv.FormatFloat(m2cm(xStep*float64(xIndex)), 'f', -1, 64), strconv.FormatFloat(electronDensity[xIndex]/(1.e6), 'f', -1, 64)})
			}
			w := csv.NewWriter(file)
			w.WriteAll(output)
			println("Electron density saved")
			if err := w.Error(); err != nil {
				log.Fatalln("error writing csv:", err)
			}
		}
	}

	if *saveSourceTermFormula {
		file, err := os.Create(configFileName + "_st.txt")
		if err != nil {
			println("unable to save ionization source term: ", err)
		} else {
			output := [][]string{{"px (Torr cm)", "a/p (cm^-1 Torr^-1)"}}
			for xIndex := range meanEnergy {
				output = append(output, []string{strconv.FormatFloat(m2cm(xStep*float64(xIndex)), 'f', -1, 64), strconv.FormatFloat(normalizedIonizationSourceTerm[xIndex]/100., 'f', -1, 64)})
			}
			w := csv.NewWriter(file)
			w.WriteAll(output)
			println("Source term saved")
			if err := w.Error(); err != nil {
				log.Fatalln("error writing csv:", err)
			}
		}
	}
	fmt.Printf("Elapsed time: %v\n", time.Since(startTime))
}
