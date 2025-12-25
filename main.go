package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"runtime"
	"runtime/pprof"
	"strings"
	"time"

	"github.com/facette/natsort"
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
	natsort.Sort(groupNames)
	for _, gName := range groupNames {
		modelNames := groups[gName]
		natsort.Sort(modelNames)

		dataRows := make(map[config.CalculationMode][]utils.Indexable[float64])
		altDataRows := make(map[config.CalculationMode][]utils.Indexable[float64])

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

			csID := parameters.CrossSections + parameters.ElasticScatteringMode + c.InelasticScatteringMode + config.StringifyMixture(parameters.GetMixtureParameters())
			if _, exists := speciesCrossSections[csID]; !exists {
				elasticSM, inelasticSM := parameters.GetScatteringModes()
				if crossSections, err := lxgata.LoadCrossSections(
					parameters.CrossSections,
					true,
					0.001,
					1000,
					elasticSM,
					inelasticSM,
					lxgata.Hartree,
					lxgata.IgnoreAtomicNumber,
					parameters.GetMixtureParameters(),
				); err == nil {
					speciesCrossSections[csID] = &crossSections
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

			var dataRow, altRow utils.Indexable[float64]
			switch parameters.CalculationMode {
			case config.BasicCalculation:
				dataRow, altRow, m, altM = extensions.BasicCalculation(parameters, c.OutputDir, modelNames[i])
			case config.CurrentCalculation:
				dataRow, altRow, m, altM = extensions.CurrentDensityCalculation(parameters, c.OutputDir, modelNames[i])
			case config.GammaCalculation:
				dataRow, altRow, m, altM = extensions.SourceIntegralCalculation(parameters, c.OutputDir, modelNames[i])
			case config.ReturningElectronsCalculation:
			case config.VoltageCalculation:
				dataRow, altRow, m, altM = extensions.VoltageCalculation2(parameters, c.OutputDir, modelNames[i])
			default:
				panic(fmt.Sprintf("unexpected config.CalculationMode: %#v", parameters.CalculationMode))
			}
			dataRow.(*utils.CoreResult).ModelName = m.Parameters.PrototypeName()
			config.AnyToSIeV(dataRow, parameters.OutputUnits(), false)
			dataRows[parameters.CalculationMode] = append(dataRows[parameters.CalculationMode], dataRow)
			if altRow != nil && altM != nil {
				altRow.(*utils.CoreResult).ModelName = altM.Parameters.PrototypeName()
				config.AnyToSIeV(altRow, parameters.OutputUnits(), false)
				altDataRows[parameters.CalculationMode] = append(altDataRows[parameters.CalculationMode], altRow)
			}

			// if parameters.CalculateSecondaryEmissionCoefficient {
			// 	// var gammaDataRow *extensions.GammaDataRow
			// 	// gammaDataRow, m = extensions.GammaCalculation(configFlags, parameters, c.OutputDir)
			// 	// gammaDataRow.ModelName = modelName
			// 	// config.AnyToSIeV(gammaDataRow, []string{"CathodeFallLength"}, parameters.OutputUnits(), false)
			// 	// gammaData = append(gammaData, gammaDataRow)
			// 	var siDataRow extensions.SourceIntegralDataRow
			// 	var altRow *extensions.SourceIntegralDataRow
			// 	siDataRow, altRow, m, altM = extensions.SourceIntegralCalculation(parameters, c.OutputDir, modelNames[i])
			// 	siDataRow.ModelName = m.Parameters.PrototypeName()
			// 	config.AnyToSIeV(&siDataRow, parameters.OutputUnits(), false)
			// 	sourceIntegralData = append(sourceIntegralData, &siDataRow)
			// 	if altRow != nil {
			// 		altRow.ModelName = altM.Parameters.PrototypeName()
			// 		config.AnyToSIeV(altRow, parameters.OutputUnits(), false)
			// 		altSourceIntegralData = append(altSourceIntegralData, altRow)
			// 	}
			// } else if parameters.CalculateCurrentDensity {
			// 	var currentDataRow extensions.CurrentDensityDataRow
			// 	var altRow *extensions.CurrentDensityDataRow
			// 	currentDataRow, altRow, m, altM = extensions.CurrentDensityCalculation(parameters, c.OutputDir, modelNames[i])
			// 	currentDataRow.ModelName = m.Parameters.PrototypeName()
			// 	config.AnyToSIeV(&currentDataRow, parameters.OutputUnits(), false)
			// 	currentData = append(currentData, &currentDataRow)
			// 	if altRow != nil {
			// 		altRow.ModelName = altM.Parameters.PrototypeName()
			// 		config.AnyToSIeV(altRow, parameters.OutputUnits(), false)
			// 		altCurrentData = append(altCurrentData, altRow)
			// 	} else {
			// 		altCurrentData = append(altCurrentData, &currentDataRow)
			// 	}
			// } else if parameters.CalculateVoltage {
			// 	var dataRow extensions.VoltageDataRow
			// 	var altRow *extensions.VoltageDataRow
			// 	dataRow, altRow, m, altM = extensions.VoltageCalculation2(parameters, c.OutputDir, modelNames[i])
			// 	dataRow.ModelName = m.Parameters.PrototypeName()
			// 	config.AnyToSIeV(&dataRow, parameters.OutputUnits(), false)
			// 	voltageData = append(voltageData, &dataRow)
			// 	if altRow != nil {
			// 		altRow.ModelName = altM.Parameters.PrototypeName()
			// 		config.AnyToSIeV(altRow, parameters.OutputUnits(), false)
			// 		altVoltageData = append(altVoltageData, altRow)
			// 	}
			// } else {
			// 	m = model.NewModel(parameters)

			// 	extensions.LoadExtensions(m.DataHub)
			// 	// additionStep := 0
			// 	m.Run(func(m *model.Model) int {
			// 		if m.TotalElectronsEmittedOnCathode == 0 {
			// 			return m.Parameters.NElectrons
			// 		} else {
			// 			collisions := m.GetMetrics(extensions.SingleElectronDetailedCollisionRateKey).(map[string]utils.GriddedInterval)
			// 			for substring := range m.Parameters.RequireCollisionRelativeMargin {
			// 				for processName, counters := range collisions {
			// 					if strings.Contains(processName, substring) {
			// 						excessLeap := utils.RelativeExcessLeap(counters.Values)
			// 						if excessLeap > 8.5 {
			// 							fmt.Printf("\nrelative excess leap for %s is %f\n", processName, excessLeap)
			// 							fmt.Printf("+ %d electrons\n", m.Parameters.AddByNElectrons)
			// 							return m.Parameters.AddByNElectrons
			// 						}
			// 					}
			// 				}
			// 			}
			// 			return 0
			// 		}
			// 	})
			// }
			if m != nil {
				output.Save(modelNames[i], m, dataExtractorFlags, c.OutputDir)
			}

			if altM != nil {
				output.Save(modelNames[i]+"_ALT", altM, dataExtractorFlags, c.OutputDir)
			}

			fmt.Printf("Elapsed time: %v\n", time.Since(modelStartTime))
		}
		for mode, list := range dataRows {
			var err error
			switch mode {
			case config.BasicCalculation:
				concreteList := make([]*utils.CoreResult, len(list))
				for i, elem := range list {
					concreteList[i] = elem.(*utils.CoreResult)
				}
				err = utils.WriteAsCSV(config.MakeHeader(concreteList[0], c.OutputUnits), concreteList, c.OutputDir+"/"+gName, "result", true)
			case config.CurrentCalculation:
				concreteList := make([]*extensions.CurrentDensityDataRow, len(list))
				for i, elem := range list {
					concreteList[i] = elem.(*extensions.CurrentDensityDataRow)
				}
				err = utils.WriteAsCSV(config.MakeHeader(concreteList[0], c.OutputUnits), concreteList, c.OutputDir+"/"+gName, "result", true)
			case config.GammaCalculation:
				concreteList := make([]*extensions.SourceIntegralDataRow, len(list))
				for i, elem := range list {
					concreteList[i] = elem.(*extensions.SourceIntegralDataRow)
				}
				err = utils.WriteAsCSV(config.MakeHeader(concreteList[0], c.OutputUnits), concreteList, c.OutputDir+"/"+gName, "result", true)
			case config.ReturningElectronsCalculation:
			case config.VoltageCalculation:
				concreteList := make([]*extensions.VoltageDataRow, len(list))
				for i, elem := range list {
					concreteList[i] = elem.(*extensions.VoltageDataRow)
				}
				err = utils.WriteAsCSV(config.MakeHeader(concreteList[0], c.OutputUnits), concreteList, c.OutputDir+"/"+gName, "result", true)
			default:
				panic(fmt.Sprintf("unexpected config.CalculationMode: %#v", mode))
			}

			if err != nil {
				fmt.Printf("unable to write data for %s, error: %v", c.OutputDir+"/"+gName, err)
			}
		}
		for mode, list := range altDataRows {
			var err error
			switch mode {
			case config.BasicCalculation:
				concreteList := make([]utils.CoreResult, len(list))
				for i, elem := range list {
					concreteList[i] = elem.(utils.CoreResult)
				}
				err = utils.WriteAsCSV(config.MakeHeader(concreteList[0], c.OutputUnits), list, c.OutputDir+"/"+gName, "result", true)
			case config.CurrentCalculation:
				concreteList := make([]extensions.CurrentDensityDataRow, len(list))
				for i, elem := range list {
					concreteList[i] = elem.(extensions.CurrentDensityDataRow)
				}
				err = utils.WriteAsCSV(config.MakeHeader(concreteList[0], c.OutputUnits), list, c.OutputDir+"/"+gName, "result", true)
			case config.GammaCalculation:
				concreteList := make([]extensions.SourceIntegralDataRow, len(list))
				for i, elem := range list {
					concreteList[i] = elem.(extensions.SourceIntegralDataRow)
				}
				err = utils.WriteAsCSV(config.MakeHeader(concreteList[0], c.OutputUnits), list, c.OutputDir+"/"+gName, "result", true)
			case config.ReturningElectronsCalculation:
			case config.VoltageCalculation:
				concreteList := make([]extensions.VoltageDataRow, len(list))
				for i, elem := range list {
					concreteList[i] = elem.(extensions.VoltageDataRow)
				}
				err = utils.WriteAsCSV(config.MakeHeader(concreteList[0], c.OutputUnits), list, c.OutputDir+"/"+gName, "result", true)
			default:
				panic(fmt.Sprintf("unexpected config.CalculationMode: %#v", mode))
			}

			if err != nil {
				fmt.Printf("unable to write data for %s, error: %v", c.OutputDir+"/"+gName, err)
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
