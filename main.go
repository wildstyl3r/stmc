package main

import (
	"flag"
	"fmt"
	"log"
	"math"
	"os"
	"runtime"
	"runtime/pprof"
	"strings"
	"time"

	"github.com/wildstyl3r/lxgata"
)

//AINF
// AIV?

func main() {
	cpuprofile := flag.String("cpuprofile", "", "write cpu profile to file")
	threads := flag.Int("j", runtime.NumCPU(), "threads to run")
	dataExtractor := newDataExtractor()
	configFileNamePointer := flag.String("input", "He_Tran_norm", "model configuration in toml format")
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

	var config, meta = loadConfig(strings.TrimSuffix(*configFileNamePointer, ".toml"))

	if config.OutputDir != "" && config.OutputDir != "." {
		os.MkdirAll(config.OutputDir, 0750)
		dataExtractor.outputPath += config.OutputDir + "/"
	}

	for modelName, parameters := range config.Models {
		fmt.Println("\n" + modelName)
		parameters.checkDefaults(modelName, &config, &meta)

		crossSections, err := lxgata.LoadCrossSections(parameters.CrossSections)
		if err != nil {
			panic(fmt.Errorf("invalid cross section file: %w", err))
		}

		model := Model{
			crossSections:     crossSections,
			cathodeFallLength: cm2m(parameters.CathodeFallLength),
			Vc:                -math.Abs(parameters.CathodeFallPotential),
			gapLength:         cm2m(parameters.GapLength),
			EConst:            parameters.ConstEField,
			density:           p / (k * parameters.Temperature),
			nElectrons:        parameters.NElectrons,
			eStep:             parameters.DeltaE,
			muStep:            parameters.DeltaMu,
			threads:           *threads,
		}
		model.init()
		model.run()
		dataExtractor.extract(modelName, &model, &parameters)
	}
	fmt.Printf("Elapsed time: %v\n", time.Since(startTime))
}
