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

func main() {
	cpuprofile := flag.String("cpuprofile", "", "write cpu profile to file")
	threads := flag.Int("j", runtime.NumCPU(), "threads to run")
	dataExtractorFlags := newDataFlags()
	configFileNamePointer := flag.String("input", "inputs/val/BM_He", "model configuration in toml format") //"inputs/val/BM_He", "model configuration in toml format")
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
		dataExtractorFlags.outputPath += config.OutputDir + "/"
	}

	for modelName, parameters := range config.Models {
		runtime.GC()
		fmt.Println("\n" + modelName)
		parameters.checkMetrizeSetDefaults(modelName, &config, &meta)

		crossSections, err := lxgata.LoadCrossSections(parameters.CrossSections)
		if err != nil {
			panic(fmt.Errorf("invalid cross section file: %w", err))
		}

		model := Model{
			crossSections: crossSections,
			dataFlags:     dataExtractorFlags,
			Vc:            -math.Abs(parameters.CathodeFallPotential),
			parameters:    parameters,
			gasDensity:    parameters.Pressure / (k * parameters.Temperature),
			threads:       *threads,
		}
		model.init()
		model.run()
		dataExtractor := DataExtractor{
			flags: dataExtractorFlags,
		}
		dataExtractor.extract(modelName, &model, &parameters)
	}
	fmt.Printf("Elapsed time: %v\n", time.Since(startTime))
}
