package main

import (
	"encoding/csv"
	"flag"
	"fmt"
	"log"
	"os"
	"sort"
	"strconv"

	"github.com/wildstyl3r/lxgata"
)

type DataItem struct {
	saveFlag    *bool
	fileSuffix  string
	columnNames []string
	values      func(*DataExtractor, int) (float64, []float64)
	scalers     []func(float64) float64
}

type DataExtractor struct {
	crossSections                    lxgata.Collisions
	xStep                            float64
	Navg, Vxavg, Favg, Eavg          float64
	cathodeFlux                      float64
	all                              *bool
	specification                    map[string]DataItem
	outputPath                       string
	TownsendAlphaF                   []float64
	rateIntegral                     []float64
	electronDensity                  []float64
	meanEnergy                       []float64
	ionizations                      []float64
	nulls                            []float64
	elastics                         []float64
	excitations                      []float64
	flux                             []float64
	rawFlux                          []float64
	energyLossByProcess              [][]float64
	probabilitiesByProcess           [][]float64
	probabilitiesByProcessFromEnergy [][]float64
}

func newDataExtractor() DataExtractor {
	return DataExtractor{
		all: flag.Bool("all", false, "save every available metric"),
		specification: map[string]DataItem{
			"Actual process probabilities": {
				saveFlag:    flag.Bool("probs", false, "save in-simulation process probabilities"),
				fileSuffix:  "probs",
				columnNames: []string{"px (Torr cm)", "..."},
				values: func(de *DataExtractor, x int) (float64, []float64) {
					return de.xStep * float64(x), de.probabilitiesByProcess[x]
				},
				scalers: []func(float64) float64{m2cm},
			},
			"Process probabilities from mean energy": {
				saveFlag:    flag.Bool("meprobs", false, "save process probabilities based on mean energy"),
				fileSuffix:  "meprobs",
				columnNames: []string{"px (Torr cm)", "..."},
				values: func(de *DataExtractor, x int) (float64, []float64) {
					return de.xStep * float64(x), de.probabilitiesByProcessFromEnergy[x]
				},
				scalers: []func(float64) float64{m2cm},
			},
			"Townsend alpha": {
				saveFlag:    flag.Bool("ta", false, "save Townsend alpha coefficient"),
				fileSuffix:  "Townsend_alpha",
				columnNames: []string{"px (Torr cm)", "a/p (cm^-1 Torr^-1)"},
				values: func(de *DataExtractor, x int) (float64, []float64) {
					return de.xStep * float64(x), []float64{de.rateIntegral[x] * de.electronDensity[x] / de.flux[x]}
				},
				scalers: []func(float64) float64{m2cm, cm2m},
			},
			"Townsend alpha per density": {
				saveFlag:    flag.Bool("tad", false, "save Townsend alpha coefficient divided by electron density"),
				fileSuffix:  "Townsend_alpha_D",
				columnNames: []string{"px (Torr cm)", "a*n/p (cm^2 Torr^-1)"},
				values: func(de *DataExtractor, x int) (float64, []float64) {
					return de.xStep * float64(x), []float64{de.rateIntegral[x] / de.flux[x]}
				},
				scalers: []func(float64) float64{m2cm, func(f float64) float64 { return f * 1.e4 }},
			},
			"Drift velocity": {
				saveFlag:    flag.Bool("vx", false, "save drift velocity"),
				fileSuffix:  "drift_vel",
				columnNames: []string{"px (Torr cm)", "v_x (cm s^-1)"},
				values: func(de *DataExtractor, x int) (float64, []float64) {
					return de.xStep * float64(x), []float64{de.flux[x] / de.electronDensity[x]}
				},
				scalers: []func(float64) float64{m2cm, m2cm},
			},
			"Density": {
				saveFlag:    flag.Bool("n", false, "save density"),
				fileSuffix:  "density",
				columnNames: []string{"px (Torr cm)", "n (сm^-3)"},
				values: func(de *DataExtractor, x int) (float64, []float64) {
					return de.xStep * float64(x), []float64{de.electronDensity[x]}
				},
				scalers: []func(float64) float64{m2cm, func(f float64) float64 { return f / 1.e6 }},
			},
			"Mean energy": {
				saveFlag:    flag.Bool("emean", false, "save mean energy"),
				fileSuffix:  "mean_energy",
				columnNames: []string{"px (Torr cm)", "e (eV)"},
				values: func(de *DataExtractor, x int) (float64, []float64) {
					return de.xStep * float64(x), []float64{de.meanEnergy[x]}
				},
				scalers: []func(float64) float64{m2cm},
			},
			"Ionization source term": {
				saveFlag:    flag.Bool("stf", false, "save ionization source term"),
				fileSuffix:  "ist",
				columnNames: []string{"px (Torr cm)", "S_n/p (cm^-1 Torr^-1)"},
				values: func(de *DataExtractor, x int) (float64, []float64) {
					return de.xStep * float64(x), []float64{de.rateIntegral[x] * de.electronDensity[x] / de.cathodeFlux}
				},
				scalers: []func(float64) float64{m2cm, cm2m},
			},
			"Ionization source term per density": {
				saveFlag:    flag.Bool("std", false, "save ionization source term divided by electron density"),
				fileSuffix:  "istd",
				columnNames: []string{"px (Torr cm)", "S_n*n/p (cm^2 Torr^-1)"},
				values: func(de *DataExtractor, x int) (float64, []float64) {
					return de.xStep * float64(x), []float64{de.rateIntegral[x] / de.cathodeFlux}
				},
				scalers: []func(float64) float64{m2cm, func(f float64) float64 { return f * 1.e4 }},
			},
			"Ionization source term per density, raw": {
				saveFlag:    flag.Bool("stdr", false, "save ionization source term with flux from naive counter, divided by electron density"),
				fileSuffix:  "istdr",
				columnNames: []string{"px (Torr cm)", "S_n*n/p (cm^2 Torr^-1)"},
				values: func(de *DataExtractor, x int) (float64, []float64) {
					return de.xStep * float64(x), []float64{de.rateIntegral[x] / (de.flux[x] / de.electronDensity[x]) * de.rawFlux[x] / de.cathodeFlux}
				},
				scalers: []func(float64) float64{m2cm, func(f float64) float64 { return f * 1.e4 }},
			},
			"Ionization counter": {
				saveFlag:    flag.Bool("ic", false, "save ionization counter"),
				fileSuffix:  "ic",
				columnNames: []string{"px (Torr cm)", "N_i(cm^-1 Torr^-1)"},
				values: func(de *DataExtractor, x int) (float64, []float64) {
					return de.xStep * float64(x), []float64{de.ionizations[x]}
				},
				scalers: []func(float64) float64{m2cm},
			},
			"Null counter": {
				saveFlag:    flag.Bool("0c", false, "save null collision counter"),
				fileSuffix:  "0c",
				columnNames: []string{"px (Torr cm)", "N_i(cm^-1 Torr^-1)"},
				values: func(de *DataExtractor, x int) (float64, []float64) {
					return de.xStep * float64(x), []float64{de.nulls[x]}
				},
				scalers: []func(float64) float64{m2cm},
			},
			"Excitations counter": {
				saveFlag:    flag.Bool("xc", false, "save excitations counter"),
				fileSuffix:  "xc",
				columnNames: []string{"px (Torr cm)", "N_i(cm^-1 Torr^-1)"},
				values: func(de *DataExtractor, x int) (float64, []float64) {
					return de.xStep * float64(x), []float64{de.excitations[x]}
				},
				scalers: []func(float64) float64{m2cm},
			},
			"Elastics counter": {
				saveFlag:    flag.Bool("ec", false, "save elastics counter"),
				fileSuffix:  "ec",
				columnNames: []string{"px (Torr cm)", "N_i(cm^-1 Torr^-1)"},
				values: func(de *DataExtractor, x int) (float64, []float64) {
					return de.xStep * float64(x), []float64{de.elastics[x]}
				},
				scalers: []func(float64) float64{m2cm},
			},
			"Electron flux": {
				saveFlag:    flag.Bool("f", false, "save calculated flux"),
				fileSuffix:  "flux",
				columnNames: []string{"px (Torr cm)", "Phi(cm^-1 Torr^-1)"},
				values: func(de *DataExtractor, x int) (float64, []float64) {
					return de.xStep * float64(x), []float64{de.flux[x]}
				},
				scalers: []func(float64) float64{m2cm, cm2m},
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
				values: func(de *DataExtractor, x int) (float64, []float64) {
					return de.xStep * float64(x), []float64{de.rateIntegral[x]}
				},
				scalers: []func(float64) float64{m2cm},
			},
			"Ionization rate mul density": {
				saveFlag:    flag.Bool("ird", false, "save ionization rate mul density"),
				fileSuffix:  "irxd",
				columnNames: []string{"px (Torr cm)", "R (cm^-3 s^-1)"},
				values: func(de *DataExtractor, x int) (float64, []float64) {
					return de.xStep * float64(x), []float64{de.rateIntegral[x] * de.electronDensity[x]}
				},
				scalers: []func(float64) float64{m2cm, func(f float64) float64 { return f / 1.e6 }},
			},
			"Ionization rate per flux": {
				saveFlag:    flag.Bool("irf", false, "save ionization rate per flux"),
				fileSuffix:  "irpf",
				columnNames: []string{"px (Torr cm)", "RpF (cm^2)"},
				values: func(de *DataExtractor, x int) (float64, []float64) {
					return de.xStep * float64(x), []float64{de.rateIntegral[x] / de.flux[x]}
				},
				scalers: []func(float64) float64{m2cm, func(f float64) float64 { return f * 1.e4 }},
			},
			"Ionization rate per velocity": {
				saveFlag:    flag.Bool("irpv", false, "save ionization rate per velocity"),
				fileSuffix:  "irpv",
				columnNames: []string{"px (Torr cm)", "RpV (cm^-1)"},
				values: func(de *DataExtractor, x int) (float64, []float64) {
					return de.xStep * float64(x), []float64{de.rateIntegral[x] * de.electronDensity[x] / de.flux[x]}
				},
				scalers: []func(float64) float64{m2cm, cm2m},
			},
			"Ionization rate per density": {
				saveFlag:    flag.Bool("irpd", false, "save ionization rate per density"),
				fileSuffix:  "irpd",
				columnNames: []string{"px (Torr cm)", "RpV (cm^3 s^-1)"},
				values: func(de *DataExtractor, x int) (float64, []float64) {
					return de.xStep * float64(x), []float64{de.rateIntegral[x] / de.electronDensity[x]}
				},
				scalers: []func(float64) float64{m2cm, func(f float64) float64 { return f * 1.e6 }},
			},
			"I*Navg": {
				saveFlag:    flag.Bool("ain", false, "i*navg"),
				fileSuffix:  "ain",
				columnNames: []string{"px (Torr cm)", "??"},
				values: func(de *DataExtractor, x int) (float64, []float64) {
					return de.xStep * float64(x), []float64{de.rateIntegral[x] * de.Navg}
				},
				scalers: []func(float64) float64{m2cm},
			},
			"I/Vavg": {
				saveFlag:    flag.Bool("aiv", false, "i/vavg"),
				fileSuffix:  "aiv",
				columnNames: []string{"px (Torr cm)", "??"},
				values: func(de *DataExtractor, x int) (float64, []float64) {
					return de.xStep * float64(x), []float64{de.rateIntegral[x] / de.Vxavg}
				},
				scalers: []func(float64) float64{m2cm},
			},
			"I/Favg": {
				saveFlag:    flag.Bool("aif", false, "i/favg"),
				fileSuffix:  "aif",
				columnNames: []string{"px (Torr cm)", "??"},
				values: func(de *DataExtractor, x int) (float64, []float64) {
					return de.xStep * float64(x), []float64{de.rateIntegral[x] / de.Favg}
				},
				scalers: []func(float64) float64{m2cm},
			},
			"I*Navg/Favg": {
				saveFlag:    flag.Bool("ainf", false, "i*navg"),
				fileSuffix:  "ainf",
				columnNames: []string{"px (Torr cm)", "??"},
				values: func(de *DataExtractor, x int) (float64, []float64) {
					return de.xStep * float64(x), []float64{de.rateIntegral[x] * de.Navg / de.Favg}
				},
				scalers: []func(float64) float64{m2cm},
			},
			"I/VavgFx": {
				saveFlag:    flag.Bool("aivx", false, "i/vavg*f(x)"),
				fileSuffix:  "aivx",
				columnNames: []string{"px (Torr cm)", "??"},
				values: func(de *DataExtractor, x int) (float64, []float64) {
					return de.xStep * float64(x), []float64{de.rateIntegral[x] / de.Vxavg * de.flux[x]}
				},
				scalers: []func(float64) float64{m2cm},
			},
			"I*Navg/FavgFx": {
				saveFlag:    flag.Bool("ainfx", false, "i*navg/favg*f(x)"),
				fileSuffix:  "ainfx",
				columnNames: []string{"px (Torr cm)", "??"},
				values: func(de *DataExtractor, x int) (float64, []float64) {
					return de.xStep * float64(x), []float64{de.rateIntegral[x] * de.Navg / de.Favg * de.flux[x]}
				},
				scalers: []func(float64) float64{m2cm},
			},
			"Energy loss due to ionizations": {
				saveFlag:    flag.Bool("li", false, "save ionization energy losses"),
				fileSuffix:  "li",
				columnNames: []string{"eV", "cm ^ -1"},
				values: func(de *DataExtractor, x int) (float64, []float64) {
					return de.xStep * float64(x), []float64{de.energyLossByProcess[x][0]}
				},
				scalers: []func(float64) float64{m2cm},
			},
			"Energy loss due to elastic collisions": {
				saveFlag:    flag.Bool("le", false, "save elastic energy losses"),
				fileSuffix:  "le",
				columnNames: []string{"eV", "cm ^ -1"},
				values: func(de *DataExtractor, x int) (float64, []float64) {
					return de.xStep * float64(x), []float64{de.energyLossByProcess[x][1]}
				},
				scalers: []func(float64) float64{m2cm},
			},
			"Energy loss due to excitations": {
				saveFlag:    flag.Bool("lx", false, "save excitation energy losses"),
				fileSuffix:  "lx",
				columnNames: []string{"eV", "cm ^ -1"},
				values: func(de *DataExtractor, x int) (float64, []float64) {
					return de.xStep * float64(x), []float64{de.energyLossByProcess[x][2]}
				},
				scalers: []func(float64) float64{m2cm},
			},
		},
	}
}

func (de *DataExtractor) extract(modelName string, model *Model, parameters *ModelParameters) {
	if de.outputPath != "" && de.outputPath[len(de.outputPath)-1] != '/' {
		de.outputPath += "/"
	}

	de.crossSections = model.crossSections

	de.cathodeFlux = parameters.CathodeCurrent / electronCharge // [m^-2 s^-1]
	de.electronDensity = make([]float64, model.numCells)
	de.rateIntegral = make([]float64, model.numCells)
	de.meanEnergy = make([]float64, model.numCells)
	de.flux = make([]float64, model.numCells)
	de.rawFlux = make([]float64, model.numCells)

	de.TownsendAlphaF = make([]float64, model.numCells)
	de.ionizations = make([]float64, model.numCells)
	de.elastics = make([]float64, model.numCells)
	de.nulls = make([]float64, model.numCells)
	de.excitations = make([]float64, model.numCells)
	de.energyLossByProcess = model.energyLossByProcess
	de.probabilitiesByProcess = make([][]float64, model.numCells)
	for i := range de.probabilitiesByProcess {
		de.probabilitiesByProcess[i] = make([]float64, len(model.probabilities[i]))
	}
	de.probabilitiesByProcessFromEnergy = make([][]float64, model.numCells)
	for i := range de.probabilitiesByProcessFromEnergy {
		de.probabilitiesByProcessFromEnergy[i] = make([]float64, len(model.probabilities[i]))
	}

	psiFIncrement := de.cathodeFlux / float64(parameters.NElectrons)

	var probKeys []string
	for processType := range model.probabilities[0] {
		probKeys = append(probKeys, processType)
	}
	sort.Strings(probKeys)
	var probKeyIndex = make(map[string]int)
	for i, key := range probKeys {
		probKeyIndex[key] = i
	}

	for xIndex := 0; xIndex < model.numCells; xIndex++ {
		for processType, probs := range model.probabilities[xIndex] {
			for i := range probs {
				de.probabilitiesByProcess[xIndex][probKeyIndex[processType]] += model.probabilities[xIndex][processType][i]
			}
			de.probabilitiesByProcess[xIndex][probKeyIndex[processType]] /= float64(len(model.probabilities[xIndex][processType]))
		}
		de.rawFlux[xIndex] = float64(model.electronsAtCell[xIndex]) / float64(model.nElectrons)
		for eIndex := 1; eIndex < model.numCellsE; eIndex++ {
			currentEnergy := model.eStep * float64(eIndex)
			fXE := 0.
			for muIndex := 0; muIndex < model.numCellsMu; muIndex++ {
				currentMu := model.muStep*float64(muIndex) - 1.

				f := 0.
				if currentMu > 0.0001 {
					f = psiFIncrement / (model.lookUpVelocity[eIndex] * currentMu) * float64(model.distribution[xIndex][eIndex][muIndex])
				}

				de.electronDensity[xIndex] += f

				fXE += f

				if eIndex > 0 {
					de.meanEnergy[xIndex] += f * currentEnergy

					de.flux[xIndex] += f * model.lookUpVelocity[eIndex] * currentMu
					de.Favg += de.flux[xIndex]

				}
			}
			cs := model.crossSections.TotalCrossSectionOfKindAt(lxgata.IONIZATION, currentEnergy)
			vel := model.lookUpVelocity[eIndex]
			de.rateIntegral[xIndex] += cs * vel * fXE
		}
		de.meanEnergy[xIndex] /= de.electronDensity[xIndex]
		de.Vxavg += de.flux[xIndex] / de.electronDensity[xIndex]
		de.Navg += de.electronDensity[xIndex]
		de.Eavg += de.meanEnergy[xIndex]

		totalCS := de.crossSections.TotalCrossSectionAt(de.meanEnergy[xIndex])
		for _, collision := range de.crossSections {
			de.probabilitiesByProcessFromEnergy[xIndex][probKeyIndex[string(collision.Type)]] += collision.CrossSectionAt(de.meanEnergy[xIndex]) / totalCS
		}
	}
	de.Vxavg /= float64(model.numCells)
	de.Navg /= float64(model.numCells)
	de.Favg /= float64(model.numCells)
	de.Eavg /= float64(model.numCells)

	de.xStep = model.xStep

	// initElectronDensity := model.cathodeFlux / math.Sqrt(2.*eV2J(4.5)/me)
	ic, elc, exc, nc := 0., 0., 0., 0.
	for xIndex := 0; xIndex < model.numCells; xIndex++ {
		de.ionizations[xIndex] = float64(model.ionizationAtCell[xIndex]) / float64(model.nElectrons) // * initElectronDensity
		de.elastics[xIndex] = float64(model.elasticAtCell[xIndex]) / float64(model.nElectrons)
		de.nulls[xIndex] = float64(model.nullAtCell[xIndex]) / float64(model.nElectrons)
		de.excitations[xIndex] = float64(model.excitationAtCell[xIndex]) / float64(model.nElectrons)
		ic += float64(model.ionizationAtCell[xIndex])
		elc += float64(model.elasticAtCell[xIndex])
		nc += float64(model.nullAtCell[xIndex])
		exc += float64(model.excitationAtCell[xIndex])
	}
	ic /= float64(model.nElectrons)
	nc /= float64(model.nElectrons)
	elc /= float64(model.nElectrons)
	exc /= float64(model.nElectrons)
	fmt.Printf("Avg collisions per electron at distance %f: elastic = %f; ionization = %f; excitation = %f; null = %f", model.cathodeFallLength, elc, ic, exc, nc)

	for xIndex := 1; xIndex+2 < model.numCells; xIndex++ {
		de.TownsendAlphaF[xIndex] = (de.flux[xIndex+1] + de.flux[xIndex+2] - de.flux[xIndex] - de.flux[xIndex-1]) / (2. * de.flux[xIndex] * de.xStep)
	}

	for name, output := range de.specification {
		if *output.saveFlag || *de.all {
			var file *os.File
			var err error
			if parameters.MakeDir && output.fileSuffix != "" && output.fileSuffix != "." {
				os.MkdirAll(de.outputPath+output.fileSuffix, 0750)
				file, err = os.Create(de.outputPath + output.fileSuffix + "/" + modelName + ".txt")
			} else {
				file, err = os.Create(de.outputPath + modelName + "_" + output.fileSuffix + ".txt")
			}
			if err != nil {
				println("unable to save "+name+": ", err)
			} else {
				rows := [][]string{output.columnNames}
				for index := 0; index < model.numCells; index++ {
					xColumnValue, yColumnValues := output.values(de, index)
					var yColumnValuesStr []string
					if len(output.scalers) > 0 && output.scalers[0] != nil {
						xColumnValue = output.scalers[0](xColumnValue)
					}
					for i := range yColumnValues {
						if len(output.scalers) > 1 && output.scalers[1] != nil {
							yColumnValuesStr = append(yColumnValuesStr, strconv.FormatFloat(output.scalers[1](yColumnValues[i]), 'f', -1, 64))
						} else {
							yColumnValuesStr = append(yColumnValuesStr, strconv.FormatFloat(yColumnValues[i], 'f', -1, 64))
						}
					}

					rows = append(rows, []string{strconv.FormatFloat(xColumnValue, 'f', -1, 64)})
					rows[len(rows)-1] = append(rows[len(rows)-1], yColumnValuesStr...)
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
