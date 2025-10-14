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
	var c, meta = config.LoadConfig(configFlags)

	speciesCrossSections := make(map[string]*lxgata.Collisions)
	var gammaData []*extensions.GammaDataRow
	var currentData []*extensions.CurrentDensityDataRow
	var voltageData []*extensions.VoltageDataRow

	for modelName, parameters := range c.Models {
		modelStartTime := time.Now()
		fmt.Printf("Current time: %s\n", modelStartTime.UTC().Format(time.UnixDate))

		runtime.GC()
		fmt.Println("\n" + modelName)
		if !parameters.CheckAndUnify(modelName, &c, &meta) {
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

		var m *model.Model

		if parameters.CalculateCathodeFallLength {
			var gammaDataRow *extensions.GammaDataRow
			gammaDataRow, m = extensions.GammaCalculation(configFlags, parameters, c.OutputDir)
			gammaDataRow.ModelName = modelName
			config.AnyToSIeV(gammaDataRow, []string{"CathodeFallLength"}, parameters.OutputUnits(), false)
			gammaData = append(gammaData, gammaDataRow)
		} else if parameters.CalculateCurrentDensity {
			var currentDataRow *extensions.CurrentDensityDataRow
			currentDataRow, m = extensions.CurrentDensityCalculation(configFlags, parameters, c.OutputDir)
			currentDataRow.ModelName = modelName
			config.AnyToSIeV(currentDataRow, []string{"CathodeFallLength", "CathodeCurrentDensity"}, parameters.OutputUnits(), false)
			currentData = append(currentData, currentDataRow)
		} else if parameters.CalculateVoltage {
			var voltageDataRow *extensions.VoltageDataRow
			voltageDataRow, m = extensions.VoltageCalculation(configFlags, parameters, c.OutputDir)
			voltageDataRow.ModelName = modelName
			// config.AnyToSIeV(voltageDataRow, []string{"Voltage"}, parameters.OutputUnits(), false)
			voltageData = append(voltageData, voltageDataRow)
		} else {
			m = model.NewModel(parameters)

			extensions.LoadExtensions(m.DataHub)
			m.Run()
		}
		output.Save(modelName, m, dataExtractorFlags, c.OutputDir)
		fmt.Printf("Elapsed time: %v\n", time.Since(modelStartTime))
	}
	if len(gammaData) > 0 {
		err := utils.WriteAsCSV(gammaData, c.OutputDir, "gamma")

		if err != nil {
			println("unable to write gamma data:", err)
		}
	}
	if len(currentData) > 0 {
		err := utils.WriteAsCSV(currentData, c.OutputDir, "current")

		if err != nil {
			println("unable to write current data:", err)
		}
	}
	if len(voltageData) > 0 {
		err := utils.WriteAsCSV(voltageData, c.OutputDir, "voltage")

		if err != nil {
			println("unable to write voltage data:", err)
		}
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
