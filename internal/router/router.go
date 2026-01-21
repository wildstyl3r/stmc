package router

import (
	"fmt"
	"runtime"
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

func RunConsole(configFlags *config.Flags, dataExtractorFlags *output.DataFlags) {
	fmt.Printf("CONFIG: %s\n", *configFlags.ConfigFileNamePointer)
	var c, meta = config.LoadConfig(*configFlags)

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

		dataRows := make(map[config.CalculationMode][]utils.ResultInterface)
		altDataRows := make(map[config.CalculationMode][]utils.ResultInterface)

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

			var dataRow, altRow utils.ResultInterface
			calculationMode := parameters.GetCalculationMode()
			switch calculationMode {
			case config.BasicCalculation:
				dataRow, altRow, m, altM = extensions.BasicCalculation(parameters, c.OutputDir, modelNames[i])
			case config.CurrentCalculation:
				dataRow, altRow, m, altM = extensions.CurrentDensityCalculation(parameters, c.OutputDir, modelNames[i])
			case config.GammaCalculation:
				dataRow, altRow, m, altM = extensions.SourceIntegralCalculation(parameters, c.OutputDir, modelNames[i])
			case config.VoltageCalculation:
				dataRow, altRow, m, altM = extensions.VoltageCalculation2(parameters, c.OutputDir, modelNames[i])
			default:
				panic(fmt.Sprintf("unexpected config.CalculationMode: %#v", parameters.CalculationMode))
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
				err = utils.WriteAsCSV(config.MakeHeader(concreteList[0], c.OutputUnits), list, c.OutputDir+"/"+gName, "result", true)
			case config.CurrentCalculation:
				concreteList := make([]*extensions.CurrentDensityDataRow, len(list))
				for i, elem := range list {
					concreteList[i] = elem.(*extensions.CurrentDensityDataRow)
				}
				err = utils.WriteAsCSV(config.MakeHeader(concreteList[0], c.OutputUnits), list, c.OutputDir+"/"+gName, "result_j", true)
			case config.GammaCalculation:
				concreteList := make([]*extensions.SourceIntegralDataRow, len(list))
				for i, elem := range list {
					concreteList[i] = elem.(*extensions.SourceIntegralDataRow)
				}
				err = utils.WriteAsCSV(config.MakeHeader(concreteList[0], c.OutputUnits), list, c.OutputDir+"/"+gName, "result_gamma", true)
			case config.VoltageCalculation:
				concreteList := make([]*extensions.VoltageDataRow, len(list))
				for i, elem := range list {
					concreteList[i] = elem.(*extensions.VoltageDataRow)
				}
				err = utils.WriteAsCSV(config.MakeHeader(concreteList[0], c.OutputUnits), list, c.OutputDir+"/"+gName, "result_V", true)
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
				concreteList := make([]*utils.CoreResult, len(list))
				for i, elem := range list {
					concreteList[i] = elem.(*utils.CoreResult)
				}
				err = utils.WriteAsCSV(config.MakeHeader(concreteList[0], c.OutputUnits), list, c.OutputDir+"/"+gName, "result", true)
			case config.CurrentCalculation:
				concreteList := make([]*extensions.CurrentDensityDataRow, len(list))
				for i, elem := range list {
					concreteList[i] = elem.(*extensions.CurrentDensityDataRow)
				}
				err = utils.WriteAsCSV(config.MakeHeader(concreteList[0], c.OutputUnits), list, c.OutputDir+"/"+gName, "result_alt_j", true)
			case config.GammaCalculation:
				concreteList := make([]*extensions.SourceIntegralDataRow, len(list))
				for i, elem := range list {
					concreteList[i] = elem.(*extensions.SourceIntegralDataRow)
				}
				err = utils.WriteAsCSV(config.MakeHeader(concreteList[0], c.OutputUnits), list, c.OutputDir+"/"+gName, "result_alt_gamma", true)
			case config.VoltageCalculation:
				concreteList := make([]*extensions.VoltageDataRow, len(list))
				for i, elem := range list {
					concreteList[i] = elem.(*extensions.VoltageDataRow)
				}
				err = utils.WriteAsCSV(config.MakeHeader(concreteList[0], c.OutputUnits), list, c.OutputDir+"/"+gName, "result_alt_V", true)
			default:
				panic(fmt.Sprintf("unexpected config.CalculationMode: %#v", mode))
			}

			if err != nil {
				fmt.Printf("unable to write data for %s, error: %v", c.OutputDir+"/"+gName, err)
			}
		}
	}
}
