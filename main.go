package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"runtime"
	"runtime/pprof"
	"time"

	"github.com/wildstyl3r/lxgata"
	"github.com/wildstyl3r/stmc/internal/config"
	"github.com/wildstyl3r/stmc/internal/extensions"
	"github.com/wildstyl3r/stmc/internal/model"
	"github.com/wildstyl3r/stmc/internal/output"
	"github.com/wildstyl3r/stmc/internal/utils"
)

func main() {
	startTime := time.Now()
	configFlags := config.NewConfigFlags()
	dataExtractorFlags := output.NewDataFlags()
	flag.Parse()

	if *configFlags.CPUprofile != "" {
		f, err := os.Create(*configFlags.CPUprofile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}
	fmt.Printf("CONFIG: %s\n", *configFlags.ConfigFileNamePointer)
	var config, meta = config.LoadConfig(configFlags)

	speciesCrossSections := make(map[string]*lxgata.Collisions)
	var gammaData utils.CSV

	for modelName, parameters := range config.Models {
		modelStartTime := time.Now()
		fmt.Printf("Current time: %s\n", modelStartTime.UTC().Format(time.UnixDate))

		runtime.GC()
		fmt.Println("\n" + modelName)
		if !parameters.CheckAndUnify(modelName, &config, &meta) {
			fmt.Printf("found a problem in the config for [%v], skipping\n", modelName)
			continue
		}

		if _, exists := speciesCrossSections[parameters.CrossSections+parameters.ScatteringMode]; !exists {
			if crossSections, err := lxgata.LoadCrossSections(parameters.CrossSections, true, parameters.GetScatteringMode(), lxgata.Hartree); err == nil {
				speciesCrossSections[parameters.CrossSections] = &crossSections
			} else {
				panic(fmt.Errorf("invalid cross section file: %w", err))
			}
		}
		parameters.SetCrossSectionsData(speciesCrossSections[parameters.CrossSections])
		if *dataExtractorFlags.Sequentials[output.MeanEnergy].DataItem.SaveFlag ||
			*dataExtractorFlags.Sequentials[output.MeanVelocityX].DataItem.SaveFlag {
			parameters.CalculateDistribution = true
		}

		if parameters.CalculateCathodeFallLength {
			gammaData = append(gammaData, extensions.GammaCalculation(configFlags, parameters, modelName, config.OutputDir))
		} else {
			m := model.NewModel(parameters)

			extensions.LoadExtensions(m.DataHub)
			m.Run()
			// dataExtractor := output.NewDataExtractor(&m)

			output.Save(modelName, &m, dataExtractorFlags, config.OutputDir)
		}
		fmt.Printf("Elapsed time: %v\n", time.Since(modelStartTime))
	}
	if len(gammaData) > 0 {
		utils.WriteAsCSV(
			gammaData,
			config.OutputDir, "gamma", *configFlags.ConfigFileNamePointer,
			[]string{"E/N", "integrated secondary emission coefficient", "analytic secondary emission coefficient", "final gamma loss", "sheath length", "integrated secondary emission coefficient_conf_interval", "last step integrated SEC", "last step SEC relative error"},
		)
	}
	fmt.Printf("Total elapsed time: %v\n\n", time.Since(startTime))
	if *configFlags.MEMprofile != "" {
		f, err := os.Create(*configFlags.MEMprofile)
		if err != nil {
			log.Fatal("could not create memory profile: ", err)
		}
		defer f.Close() // error handling omitted for example
		runtime.GC()    // get up-to-date statistics
		// Lookup("allocs") creates a profile similar to go test -memprofile.
		// Alternatively, use Lookup("heap") for a profile
		// that has inuse_space as the default index.
		if err := pprof.Lookup("allocs").WriteTo(f, 0); err != nil {
			log.Fatal("could not write memory profile: ", err)
		}
	}
}
