package router

import (
	"fmt"
	"maps"
	"runtime"
	"slices"
	"strings"
	"time"

	"github.com/facette/natsort"
	"github.com/wildstyl3r/lxgata"
	"github.com/wildstyl3r/stmc/internal/config"
	"github.com/wildstyl3r/stmc/internal/extensions"
	"github.com/wildstyl3r/stmc/internal/messages"
	"github.com/wildstyl3r/stmc/internal/model"
	"github.com/wildstyl3r/stmc/internal/output"
	"github.com/wildstyl3r/stmc/internal/utils"
)

func RunConsole(configFlags *config.Flags, dataExtractorFlags *output.DataFlags, logger messages.Logger) {
	fmt.Printf("CONFIG: %s\n", *configFlags.ConfigFileNamePointer)
	var c, meta = config.LoadConfig(*configFlags, logger)

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
	modelStatus := make(map[string]string)
	for _, gName := range groupNames {
		modelNames := groups[gName]
		natsort.Sort(modelNames)

		dataRows := make(map[config.CalculationMode][]utils.ResultInterface)
		altDataRows := make(map[config.CalculationMode][]utils.ResultInterface)

		for i := range modelNames {
			parameters := c.Models[modelNames[i]]
			modelStartTime := time.Now()
			fmt.Printf("Current time: %s\n", modelStartTime.UTC().Format(time.UnixDate))

			runtime.GC()
			fmt.Println("\n" + modelNames[i])
			if !parameters.CheckAndUnify(modelNames[i], &c, &meta, logger) {
				// fmt.Printf("found a problem in the config for [%v], skipping\n", modelNames[i])
				modelStatus[modelNames[i]] = fmt.Sprintf("found a problem in the config for [%v], skipping\n", modelNames[i])
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
					"electron",
				); err == nil {
					speciesCrossSections[csID] = &crossSections
				} else {
					modelStatus[modelNames[i]] = fmt.Sprintf("invalid cross section file: %v\n", err)
					continue
				}
			}
			parameters.SetCrossSectionsData(speciesCrossSections[csID])
			if *dataExtractorFlags.Sequentials[output.MeanEnergy].DataItem.SaveFlag ||
				*dataExtractorFlags.Sequentials[output.MeanVelocityX].DataItem.SaveFlag {
				parameters.CalculateDistribution = true
			}

			var m, altM *model.Model

			var dataRow, altRow utils.ResultInterface
			calculationMode := parameters.GetCalculationMode()
			switch calculationMode {
			case config.BasicSheathCalculation:
				dataRow, altRow, m, altM = extensions.BasicCalculation(parameters, c.OutputDir, modelNames[i], logger)
			case config.CurrentCalculation:
				dataRow, altRow, m, altM = extensions.CurrentDensityCalculation(parameters, c.OutputDir, modelNames[i], logger)
			case config.GammaCalculation:
				dataRow, altRow, m, altM = extensions.SourceIntegralCalculation(parameters, c.OutputDir, modelNames[i], logger)
			case config.VoltageCalculation:
				dataRow, altRow, m, altM = extensions.VoltageCalculation2(parameters, c.OutputDir, modelNames[i], logger)
			case config.TownsendAlpha:
				dataRow, altRow, m, altM = extensions.AlphaCalculationR(parameters, c.OutputDir, modelNames[i], logger)
			default:
				logger.Failure("unexpected config.CalculationMode: %#v", parameters.CalculationMode)
			}
			// dataRow.(*utils.CoreResult).ModelName = m.Parameters.PrototypeName()
			config.AnyToSIeV(dataRow, parameters.OutputUnits(), false)
			dataRows[calculationMode] = append(dataRows[calculationMode], dataRow)
			if altRow != nil && altM != nil {
				// altRow.(*utils.CoreResult).ModelName = altM.Parameters.PrototypeName()
				config.AnyToSIeV(altRow, parameters.OutputUnits(), false)
				altDataRows[calculationMode] = append(altDataRows[calculationMode], altRow)
			}

			if m != nil {
				output.Save(modelNames[i], m, *dataExtractorFlags, c.OutputDir)
			}

			if altM != nil {
				output.Save(modelNames[i]+"_ALT", altM, *dataExtractorFlags, c.OutputDir)
			}

			modelStatus[modelNames[i]] = "success\n"

			fmt.Printf("Elapsed time: %v\n", time.Since(modelStartTime))
		}
		for mode, list := range dataRows {
			var err error
			switch mode {
			case config.BasicSheathCalculation:
				specificList := make([]*utils.SheathResult, len(list))
				for i, elem := range list {
					specificList[i] = elem.(*utils.SheathResult)
				}
				err = utils.WriteAsCSV(config.MakeHeader(specificList[0], c.OutputUnits), list, c.OutputDir+"/"+gName, "result", true)
			case config.CurrentCalculation:
				specificList := make([]*extensions.CurrentDensityDataRow, len(list))
				for i, elem := range list {
					specificList[i] = elem.(*extensions.CurrentDensityDataRow)
				}
				err = utils.WriteAsCSV(config.MakeHeader(specificList[0], c.OutputUnits), list, c.OutputDir+"/"+gName, "result_j", true)
			case config.GammaCalculation:
				specificList := make([]*extensions.SourceIntegralDataRow, len(list))
				for i, elem := range list {
					specificList[i] = elem.(*extensions.SourceIntegralDataRow)
				}
				err = utils.WriteAsCSV(config.MakeHeader(specificList[0], c.OutputUnits), list, c.OutputDir+"/"+gName, "result_gamma", true)
			case config.VoltageCalculation:
				specificList := make([]*extensions.VoltageDataRow, len(list))
				for i, elem := range list {
					specificList[i] = elem.(*extensions.VoltageDataRow)
				}
				err = utils.WriteAsCSV(config.MakeHeader(specificList[0], c.OutputUnits), list, c.OutputDir+"/"+gName, "result_V", true)
			case config.TownsendAlpha:
				specificList := make([]*extensions.AlphaResult, len(list))
				for i, elem := range list {
					specificList[i] = elem.(*extensions.AlphaResult)
				}
				err = utils.WriteAsCSV(config.MakeHeader(specificList[0], c.OutputUnits), list, c.OutputDir+"/"+gName, "result_alpha", true)
			default:
				logger.Failure("unexpected config.CalculationMode: %#v", mode)
			}

			if err != nil {
				logger.Failure("unable to write data for %s, error: %v", c.OutputDir+"/"+gName, err)
			}
		}
		for mode, list := range altDataRows {
			var err error
			fail := false
			switch mode {
			case config.BasicSheathCalculation:
				specificList := make([]*utils.SheathResult, len(list))
				for i, elem := range list {
					specificList[i] = elem.(*utils.SheathResult)
				}
				err = utils.WriteAsCSV(config.MakeHeader(specificList[0], c.OutputUnits), list, c.OutputDir+"/"+gName, "result", true)
			case config.CurrentCalculation:
				specificList := make([]*extensions.CurrentDensityDataRow, len(list))
				for i, elem := range list {
					specificList[i] = elem.(*extensions.CurrentDensityDataRow)
				}
				err = utils.WriteAsCSV(config.MakeHeader(specificList[0], c.OutputUnits), list, c.OutputDir+"/"+gName, "result_alt_j", true)
			case config.GammaCalculation:
				specificList := make([]*extensions.SourceIntegralDataRow, len(list))
				for i, elem := range list {
					specificList[i] = elem.(*extensions.SourceIntegralDataRow)
				}
				err = utils.WriteAsCSV(config.MakeHeader(specificList[0], c.OutputUnits), list, c.OutputDir+"/"+gName, "result_alt_gamma", true)
			case config.VoltageCalculation:
				specificList := make([]*extensions.VoltageDataRow, len(list))
				for i, elem := range list {
					specificList[i] = elem.(*extensions.VoltageDataRow)
				}
				err = utils.WriteAsCSV(config.MakeHeader(specificList[0], c.OutputUnits), list, c.OutputDir+"/"+gName, "result_alt_V", true)
			case config.TownsendAlpha:
				specificList := make([]*extensions.AlphaResult, len(list))
				for i, elem := range list {
					specificList[i] = elem.(*extensions.AlphaResult)
				}
				err = utils.WriteAsCSV(config.MakeHeader(specificList[0], c.OutputUnits), list, c.OutputDir+"/"+gName, "result_altpha", true)
			default:
				fail = true
				logger.Failure("unexpected config.CalculationMode: %#v", mode)
			}

			if err != nil {
				logger.Info("unable to write data for %s, error: %v", c.OutputDir+"/"+gName, err)
			} else if !fail {
				logger.Info("Calculation completed")
			}
		}
	}

	allModelNames := slices.Sorted(maps.Keys(modelStatus))
	status := ""
	for _, mn := range allModelNames {
		status += mn + ": " + modelStatus[mn] + "\n"
	}
	if status != "" {
		logger.Info(status)
	}
}
