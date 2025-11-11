package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"runtime"
	"runtime/pprof"
	"slices"
	"strings"
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

	groups := make(map[string][]string)
	groupNames := make([]string, 0)
	for name := range c.Models {
		groupName, _ := strings.CutSuffix(name, c.Models[name].PrototypeName())
		if _, exist := groups[groupName]; !exist {
			groups[groupName] = make([]string, 0)
		}
		groups[groupName] = append(groups[groupName], name)

	}
	for gname := range groups {
		groupNames = append(groupNames, gname)
	}
	slices.Sort(groupNames)
	for _, gName := range groupNames {
		modelNames := groups[gName]
		slices.Sort(modelNames)

		var currentData []*extensions.CurrentDensityDataRow
		var voltageData []*extensions.VoltageDataRow
		var sourceIntegralData []*extensions.SourceIntegralDataRow
		var altCurrentData []*extensions.CurrentDensityDataRow
		var altVoltageData []*extensions.VoltageDataRow
		var altSourceIntegralData []*extensions.SourceIntegralDataRow

		for i := range modelNames {
			parameters := c.Models[modelNames[i]]
			modelStartTime := time.Now()
			fmt.Printf("Current time: %s\n", modelStartTime.UTC().Format(time.UnixDate))

			runtime.GC()
			fmt.Println("\n" + modelNames[i])
			if !parameters.CheckAndUnify(modelNames[i], &c, &meta) {
				fmt.Printf("found a problem in the config for [%v], skipping\n", modelNames[i])
				continue
			}

			csID := parameters.CrossSections + parameters.ElasticScatteringMode + c.InelasticScatteringMode
			if _, exists := speciesCrossSections[csID]; !exists {
				elasticSM, inelasticSM := parameters.GetScatteringModes()
				if crossSections, err := lxgata.LoadCrossSections(parameters.CrossSections, true, 0.001, 1000, elasticSM, inelasticSM, lxgata.Hartree, lxgata.IgnoreAtomicNumber); err == nil {
					speciesCrossSections[csID] = &crossSections
					for i := range crossSections.Processes {
						crossSections.Processes[i].Expand(0.01, false)
					}
				} else {
					panic(fmt.Errorf("invalid cross section file: %w", err))
				}
			}
			parameters.SetCrossSectionsData(speciesCrossSections[csID])
			if *dataExtractorFlags.Sequentials[output.MeanEnergy].DataItem.SaveFlag ||
				*dataExtractorFlags.Sequentials[output.MeanVelocityX].DataItem.SaveFlag {
				parameters.CalculateDistribution = true
			}

			var m, altM *model.Model

			if parameters.CalculateSecondaryEmissionCoefficient {
				// var gammaDataRow *extensions.GammaDataRow
				// gammaDataRow, m = extensions.GammaCalculation(configFlags, parameters, c.OutputDir)
				// gammaDataRow.ModelName = modelName
				// config.AnyToSIeV(gammaDataRow, []string{"CathodeFallLength"}, parameters.OutputUnits(), false)
				// gammaData = append(gammaData, gammaDataRow)
				var siDataRow extensions.SourceIntegralDataRow
				var altRow *extensions.SourceIntegralDataRow
				siDataRow, altRow, m, altM = extensions.SourceIntegralCalculation(parameters, c.OutputDir, modelNames[i])
				siDataRow.ModelName = modelNames[i]
				config.AnyToSIeV(&siDataRow, parameters.OutputUnits(), false)
				sourceIntegralData = append(sourceIntegralData, &siDataRow)
				if altRow != nil {
					altRow.ModelName = modelNames[i]
					config.AnyToSIeV(altRow, parameters.OutputUnits(), false)
					altSourceIntegralData = append(altSourceIntegralData, altRow)
				}
			} else if parameters.CalculateCurrentDensity {
				var currentDataRow extensions.CurrentDensityDataRow
				var altRow *extensions.CurrentDensityDataRow
				currentDataRow, altRow, m, altM = extensions.CurrentDensityCalculation(parameters, c.OutputDir, modelNames[i])
				currentDataRow.ModelName = modelNames[i]
				config.AnyToSIeV(&currentDataRow, parameters.OutputUnits(), false)
				currentData = append(currentData, &currentDataRow)
				if altRow != nil {
					altRow.ModelName = modelNames[i]
					config.AnyToSIeV(altRow, parameters.OutputUnits(), false)
					altCurrentData = append(altCurrentData, altRow)
				} else {
					altCurrentData = append(altCurrentData, &currentDataRow)
				}
			} else if parameters.CalculateVoltage {
				var voltageDataRow extensions.VoltageDataRow
				var altRow *extensions.VoltageDataRow
				voltageDataRow, altRow, m, altM = extensions.VoltageCalculation2(parameters, c.OutputDir, modelNames[i])
				voltageDataRow.ModelName = modelNames[i]
				config.AnyToSIeV(&voltageDataRow, parameters.OutputUnits(), false)
				voltageData = append(voltageData, &voltageDataRow)
				if altRow != nil {
					altRow.ModelName = modelNames[i]
					config.AnyToSIeV(altRow, parameters.OutputUnits(), false)
					altVoltageData = append(altVoltageData, altRow)
				}
			} else {
				m = model.NewModel(parameters)

				extensions.LoadExtensions(m.DataHub)
				// additionStep := 0
				m.Run(func(m *model.Model) int {
					if m.TotalElectronsEmittedOnCathode == 0 {
						return m.Parameters.NElectrons
					} else {
						// if additionStep < parameters.AdditionSteps {
						// additionStep += 1
						collisions := m.GetMetrics(extensions.SingleElectronDetailedCollisionRateKey).(map[string]utils.GriddedInterval)
						for substring := range m.Parameters.RequireCollisionRelativeMargin {
							for processName, counters := range collisions {
								if strings.Contains(processName, substring) {
									excessLeap := utils.RelativeExcessLeap(counters.Values)
									if excessLeap > 8.5 {
										fmt.Printf("\nrelative excess leap for %s is %f\n", processName, excessLeap)
										fmt.Printf("+ %d electrons\n", m.Parameters.AddByNElectrons)
										return m.Parameters.AddByNElectrons
									}
								}
							}
						}
						// }
						// if len(m.Parameters.RequireCollisionRelativeMargin) > 0 {
						// 	// lowerMargins := m.GetMetrics(extensions.SingleElectronDetailedCollisionRateLowerMarginKey).(map[string]utils.GriddedInterval)
						// 	upperMargins := m.GetMetrics(extensions.SingleElectronDetailedCollisionRateUpperMarginKey).(map[string]utils.GriddedInterval)
						// 	values := m.GetMetrics(extensions.SingleElectronDetailedCollisionRateKey).(map[string]utils.GriddedInterval)

						// }
						return 0
					}
				})
			}
			if m != nil {
				output.Save(modelNames[i], m, dataExtractorFlags, c.OutputDir)
			}

			if altM != nil {
				output.Save(modelNames[i]+"_ALT", altM, dataExtractorFlags, c.OutputDir)
			}

			fmt.Printf("Elapsed time: %v\n", time.Since(modelStartTime))
		}
		if len(currentData) > 0 {
			err := utils.WriteAsCSV(config.MakeHeader(currentData[0], c.OutputUnits), currentData, c.OutputDir+"/"+gName, "current", true)

			if err != nil {
				println("unable to write current data:", err)
			}
		}
		if len(voltageData) > 0 {
			err := utils.WriteAsCSV(config.MakeHeader(voltageData[0], c.OutputUnits), voltageData, c.OutputDir+"/"+gName, "voltage", true)

			if err != nil {
				println("unable to write voltage data:", err)
			}
		}
		if len(sourceIntegralData) > 0 {
			err := utils.WriteAsCSV(config.MakeHeader(sourceIntegralData[0], c.OutputUnits), sourceIntegralData, c.OutputDir+"/"+gName, "sourceInt", true)

			if err != nil {
				println("unable to write source integral data:", err)
			}
		}
		if len(altCurrentData) > 0 {
			err := utils.WriteAsCSV(config.MakeHeader(altCurrentData[0], c.OutputUnits), altCurrentData, c.OutputDir+"/"+gName, "alt_current", true)

			if err != nil {
				println("unable to write current data:", err)
			}
		}
		if len(altVoltageData) > 0 {
			err := utils.WriteAsCSV(config.MakeHeader(altVoltageData[0], c.OutputUnits), altVoltageData, c.OutputDir+"/"+gName, "alt_voltage", true)

			if err != nil {
				println("unable to write voltage data:", err)
			}
		}
		if len(altSourceIntegralData) > 0 {
			err := utils.WriteAsCSV(config.MakeHeader(altSourceIntegralData[0], c.OutputUnits), altSourceIntegralData, c.OutputDir+"/"+gName, "alt_sourceInt", true)

			if err != nil {
				println("unable to write source integral data:", err)
			}
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
