package main

import (
	"encoding/csv"
	"flag"
	"fmt"
	"log"
	"math"
	"os"
	"sort"
	"strconv"

	"github.com/wildstyl3r/lxgata"
)

type DataItem struct {
	saveFlag   *bool
	fileSuffix string
}

type SequentialDataItem struct {
	DataItem
	columnNames []string
	values      func(*DataExtractor) (args []float64, values [][]float64, labels []string)
	scalers     [][]func(float64) float64
}

type ScalarDataItem struct {
	DataItem
	value func(*DataExtractor) (value float64)
}

type DataFlags struct {
	all         *bool
	sequentials map[string]SequentialDataItem
	scalars     map[string]ScalarDataItem
	outputPath  string
}

type DataExtractor struct {
	crossSections                    *lxgata.Collisions
	model                            *Model
	psiFIncrement                    float64
	gasDensity                       float64
	cathodeFlux                      float64
	flags                            DataFlags
	TownsendAlphaF                   []float64
	rateIntegral                     []float64
	sourceTerm                       []float64
	electronDensity                  []float64
	meanEnergy                       []float64
	driftVelocity                    []float64
	angularDistribution              [][]float64
	energyDistribution               [][]float64
	collisions                       map[string][]float64
	flux                             []float64
	velocity                         []float64
	radialVelocity                   []float64
	energyLossByProcess              map[string][]float64
	probabilitiesByProcess           map[string][]float64
	probabilitiesByProcessFromEnergy map[string][]float64
}

func newDataFlags() DataFlags {
	return DataFlags{
		all: flag.Bool("all", true, "save every available metric"),
		sequentials: map[string]SequentialDataItem{

			// "Angular distribution": {
			// 	saveFlag:    flag.Bool("ang", false, "save potential"),
			// 	fileSuffix:  "angles",
			// 	columnNames: []string{"x (cm)", "cos(Theta)"},
			// 	values: func(de *DataExtractor) (args []float64, values [][]float64, labels []string) {
			// 		for x := 0; x < de.model.numCells; x++ {
			// 			args = append(args, de.model.xStep*float64(x))
			// 			var muRow = make([]float64, de.model.numCellsMu)
			// 			for muIndex := 0; muIndex < de.model.numCellsMu; muIndex++ {
			// 				currentMu := de.model.muStep*float64(muIndex) - 1.
			// 				accum := 0.
			// 				for eIndex := 0; eIndex < de.model.numCellsE; eIndex++ {
			// 					// currentEnergy := de.model.eStep * float64(eIndex)
			// 					accum += de.psiFIncrement / (de.model.lookUpVelocity[eIndex] * currentMu) * float64(de.model.distribution[x][eIndex][muIndex])
			// 					// de.model.distribution[x][e][mu]
			// 				}
			// 				muRow = append(muRow, accum)
			// 			}
			// 			values = append(values, muRow)
			// 		}
			// 		return args, values, nil
			// 	},
			// 	scalers: []func(float64) float64{m2cm},
			// },
			// "Energy distribution": {
			// 	saveFlag:    flag.Bool("engs", false, "save potential"),
			// 	fileSuffix:  "energies",
			// 	columnNames: []string{"x (cm)", "e, eV"},
			// 	values: func(de *DataExtractor) (args []float64, values [][]float64, labels []string) {
			// 		for x := 0; x < de.model.numCells; x++ {
			// 			args = append(args, de.model.xStep*float64(x))
			// 			var eRow = make([]float64, de.model.numCellsMu)
			// 			for muIndex := 0; muIndex < de.model.numCellsMu; muIndex++ {
			// 				currentMu := de.model.muStep*float64(muIndex) - 1.
			// 				accum := 0.
			// 				for eIndex := 0; eIndex < de.model.numCellsE; eIndex++ {
			// 					// currentEnergy := de.model.eStep * float64(eIndex)
			// 					accum += de.psiFIncrement / (de.model.lookUpVelocity[eIndex] * currentMu) * float64(de.model.distribution[x][eIndex][muIndex])
			// 					// de.model.distribution[x][e][mu]
			// 				}
			// 				muRow = append(muRow, accum)
			// 			}
			// 			values = append(values, []float64{de.model.VfromL(de.model.xStep * float64(x))})
			// 		}
			// 		return args, values, nil
			// 	},
			// 	scalers: []func(float64) float64{m2cm},
			// },
			"Potential": {
				DataItem: DataItem{
					saveFlag:   flag.Bool("v", false, "save potential"),
					fileSuffix: "V",
				},
				columnNames: []string{"x (cm)", "g (V)"},
				values: func(de *DataExtractor) (args []float64, values [][]float64, labels []string) {
					for x := 0; x < de.model.numCells; x++ {
						args = append(args, de.model.xStep*float64(x))
						values = append(values, []float64{de.model.VfromL(de.model.xStep * float64(x))})
					}
					return args, values, nil
				},
				scalers: [][]func(float64) float64{{m2cm}},
			},
			"Electric Field": {
				DataItem: DataItem{
					saveFlag:   flag.Bool("ef", false, "save Electric field"),
					fileSuffix: "Efield",
				},
				columnNames: []string{"x (cm)", "E (V/m)"},
				values: func(de *DataExtractor) (args []float64, values [][]float64, labels []string) {
					for x := 0; x < de.model.numCells; x++ {
						args = append(args, de.model.xStep*float64(x))
						values = append(values, []float64{de.model.EFieldFromL(de.model.xStep * float64(x))})
					}
					return args, values, nil
				},
				scalers: [][]func(float64) float64{{m2cm}, {cm2m}},
			},
			"LfromV": {
				DataItem: DataItem{
					saveFlag:   flag.Bool("lv", false, "save x from v"),
					fileSuffix: "lv",
				},
				columnNames: []string{"g (V)", "x (cm)"},
				values: func(de *DataExtractor) (args []float64, values [][]float64, labels []string) {
					for i := 0; i < 100000; i++ {
						v := -(de.model.Va + de.model.Vc) * float64(100000-i) / 100000.
						args = append(args, v)
						values = append(values, []float64{de.model.LfromV(v)})
					}
					return args, values, nil
				},
				scalers: [][]func(float64) float64{{}, {m2cm}},
			},
			"Actual process probabilities": {
				DataItem: DataItem{
					saveFlag:   flag.Bool("probs", false, "save in-simulation process probabilities"),
					fileSuffix: "probs",
				},
				columnNames: []string{"x (cm)", "..."},
				values: func(de *DataExtractor) (args []float64, values [][]float64, labels []string) {
					for label := range de.probabilitiesByProcess {
						labels = append(labels, label)
					}
					sort.Strings(labels)
					for x := 0; x < de.model.numCells; x++ {
						args = append(args, de.model.xStep*float64(x))
						var row []float64
						for _, label := range labels {
							row = append(row, de.probabilitiesByProcess[label][x])
						}
						values = append(values, row)
					}
					return args, values, labels
				},
				scalers: [][]func(float64) float64{{m2cm}},
			},
			"Process probabilities from mean energy": {
				DataItem: DataItem{
					saveFlag:   flag.Bool("meprobs", false, "save process probabilities based on mean energy"),
					fileSuffix: "meprobs",
				},
				columnNames: []string{"x (cm)", "..."},
				values: func(de *DataExtractor) (args []float64, values [][]float64, labels []string) {
					for label := range de.probabilitiesByProcessFromEnergy {
						labels = append(labels, label)
					}
					sort.Strings(labels)
					for x := 0; x < de.model.numCells; x++ {
						args = append(args, de.model.xStep*float64(x))
						var row []float64
						for _, label := range labels {
							row = append(row, de.probabilitiesByProcessFromEnergy[label][x])
						}
						values = append(values, row)
					}
					return args, values, labels
				},
				scalers: [][]func(float64) float64{{m2cm}},
			},
			"Townsend alpha": {
				DataItem: DataItem{
					saveFlag:   flag.Bool("ta", false, "save Townsend alpha coefficient"),
					fileSuffix: "Townsend_alpha",
				},
				columnNames: []string{"x (cm)", "a/p (cm^-1 Torr^-1)"},
				values: func(de *DataExtractor) (args []float64, values [][]float64, labels []string) {
					for x := range de.rateIntegral {
						args = append(args, de.model.xStep*float64(x))
						values = append(values, []float64{de.rateIntegral[x] / de.driftVelocity[x] * de.gasDensity})
					}
					return args, values, nil
				},
				scalers: [][]func(float64) float64{{m2cm}, {cm2m}},
			},
			"Townsend alpha alt": {
				DataItem: DataItem{
					saveFlag:   flag.Bool("taa", false, "save Townsend alpha coefficient as in Boeuf & Marode 1982"),
					fileSuffix: "Townsend_alpha_alt",
				},
				columnNames: []string{"x (cm)", "a/p (cm^2 Torr^-1)"},
				values: func(de *DataExtractor) (args []float64, values [][]float64, labels []string) {
					avgDVel := 0.
					for x := range de.driftVelocity {
						avgDVel += de.driftVelocity[x]
					}
					avgDVel /= float64(len(de.driftVelocity))
					for x := range de.rateIntegral {
						args = append(args, de.model.xStep*float64(x))
						values = append(values, []float64{de.rateIntegral[x] * de.electronDensity[x] / avgDVel})
					}
					return args, values, nil
				},
				scalers: [][]func(float64) float64{{m2cm}, {cm2m}},
			},
			"Drift velocity": {
				DataItem: DataItem{
					saveFlag:   flag.Bool("vx", false, "save drift velocity"),
					fileSuffix: "drift_vel",
				},
				columnNames: []string{"x (cm)", "v_x (cm s^-1)"},
				values: func(de *DataExtractor) (args []float64, values [][]float64, labels []string) {
					for x := range de.flux {
						args = append(args, de.model.xStep*float64(x))
						values = append(values, []float64{de.flux[x] / de.electronDensity[x]})
					}
					return args, values, nil
				},
				scalers: [][]func(float64) float64{{m2cm}, {m2cm}},
			},
			"Radial velocity": {
				DataItem: DataItem{
					saveFlag:   flag.Bool("vs", false, "save radial velocity"),
					fileSuffix: "vel_star",
				},
				columnNames: []string{"x (cm)", "v* (cm s^-1)"},
				values: func(de *DataExtractor) (args []float64, values [][]float64, labels []string) {
					for x := range de.flux {
						args = append(args, de.model.xStep*float64(x))
						values = append(values, []float64{de.radialVelocity[x]})
					}
					return args, values, nil
				},
				scalers: [][]func(float64) float64{{m2cm}, {m2cm}},
			},
			"Total velocity": {
				DataItem: DataItem{
					saveFlag:   flag.Bool("vt", false, "save total velocity"),
					fileSuffix: "vel",
				},
				columnNames: []string{"x (cm)", "v (cm s^-1)"},
				values: func(de *DataExtractor) (args []float64, values [][]float64, labels []string) {
					for x := range de.flux {
						args = append(args, de.model.xStep*float64(x))
						values = append(values, []float64{de.velocity[x]})
					}
					return args, values, nil
				},
				scalers: [][]func(float64) float64{{m2cm}, {m2cm}},
			},
			"M base": {
				DataItem: DataItem{
					saveFlag:   flag.Bool("m", false, "save total cross section dynamics"),
					fileSuffix: "m",
				},
				columnNames: []string{"x (cm)", "arb"},
				values: func(de *DataExtractor) (args []float64, values [][]float64, labels []string) {
					f := func(x float64, totEnergy float64) float64 {
						eKin := totEnergy + de.model.VfromL(x)
						potential := -(totEnergy - eKin)
						return de.model.gasDensity * de.model.crossSections.TotalCrossSectionAt(eKin) * math.Sqrt(eKin) / math.Abs(de.model.EFieldFromL(de.model.LfromV(potential)))
					}
					for i := 0; i < de.model.numCells*2; i++ {
						x := de.model.xStep * float64(i) / 2.
						args = append(args, x)
						values = append(values, []float64{
							f(x, math.Abs(de.model.parameters.CathodeFallPotential)+4.5),
							f(x, (math.Abs(de.model.parameters.CathodeFallPotential)+4.5)*0.75),
							f(x, (math.Abs(de.model.parameters.CathodeFallPotential)+4.5)*0.5),
							f(x, (math.Abs(de.model.parameters.CathodeFallPotential)+4.5)*0.25)})
					}
					return args, values, []string{"total energy: 1", "total energy: 0.75", "total energy: 0.5", "total energy: 0.25"}
				},
				scalers: [][]func(float64) float64{{m2cm}, {m2cm}},
			},
			"Density": {
				DataItem: DataItem{
					saveFlag:   flag.Bool("n", false, "save density"),
					fileSuffix: "density",
				},
				columnNames: []string{"x (cm)", "n (cm^-3)"},
				values: func(de *DataExtractor) (args []float64, values [][]float64, labels []string) {
					for x := range de.electronDensity {
						args = append(args, de.model.xStep*float64(x))
						values = append(values, []float64{de.electronDensity[x]})
					}
					return args, values, nil
				},
				scalers: [][]func(float64) float64{{m2cm}, {cm2m, cm2m, cm2m}},
			},
			"Mean energy": {
				DataItem: DataItem{
					saveFlag:   flag.Bool("emean", false, "save mean energy"),
					fileSuffix: "mean_energy",
				},
				columnNames: []string{"x (cm)", "e (eV)"},
				values: func(de *DataExtractor) (args []float64, values [][]float64, labels []string) {
					for x := range de.meanEnergy {
						args = append(args, de.model.xStep*float64(x))
						values = append(values, []float64{de.meanEnergy[x]})
					}
					return args, values, nil
				},
				scalers: [][]func(float64) float64{{m2cm}},
			},
			"Source term": {
				DataItem: DataItem{
					saveFlag:   flag.Bool("st", false, "save source term"),
					fileSuffix: "st",
				},
				columnNames: []string{"x (cm)", "S (cm^-3 s^-1)"},
				values: func(de *DataExtractor) (args []float64, values [][]float64, labels []string) {
					for x := range de.sourceTerm {
						args = append(args, de.model.xStep*float64(x))
						values = append(values, []float64{de.sourceTerm[x]})
					}
					return args, values, nil
				},
				scalers: [][]func(float64) float64{{m2cm}, {cm2m, cm2m, cm2m}},
			},
			"Normalized source term": {
				DataItem: DataItem{
					saveFlag:   flag.Bool("nst", false, "save normalized source term"),
					fileSuffix: "nst",
				},
				columnNames: []string{"px (cm Torr)", "S/p (cm^-1 Torr ^-1)"},
				values: func(de *DataExtractor) (args []float64, values [][]float64, labels []string) {
					for x := range de.sourceTerm {
						args = append(args, de.model.xStep*float64(x))
						values = append(values, []float64{de.sourceTerm[x] * Torr / de.model.parameters.Pressure / de.cathodeFlux})
					}
					return args, values, nil
				},
				scalers: [][]func(float64) float64{{m2cm}, {cm2m}},
			},
			"Collision counters": {
				DataItem: DataItem{
					saveFlag:   flag.Bool("cc", false, "save collision counters"),
					fileSuffix: "cc",
				},
				columnNames: []string{"x (cm)", "N_i(cm^-1 Torr^-1)"},
				values: func(de *DataExtractor) (args []float64, values [][]float64, labels []string) {
					for label := range de.collisions {
						labels = append(labels, label)
					}
					sort.Strings(labels)
					for x := 0; x < de.model.numCells; x++ {
						args = append(args, de.model.xStep*float64(x))
						var row []float64
						for _, label := range labels {
							row = append(row, de.collisions[label][x])
						}
						values = append(values, row)
					}
					return args, values, labels
				},
				scalers: [][]func(float64) float64{{m2cm}, {cm2m}},
			},
			"Electron flux": {
				DataItem: DataItem{
					saveFlag:   flag.Bool("f", false, "save calculated flux"),
					fileSuffix: "flux",
				},
				columnNames: []string{"x (cm)", "Phi(cm^-1 Torr^-1)"},
				values: func(de *DataExtractor) (args []float64, values [][]float64, labels []string) {
					for x := range de.flux {
						args = append(args, de.model.xStep*float64(x))
						values = append(values, []float64{de.flux[x]})
					}
					return args, values, nil
				},
				scalers: [][]func(float64) float64{{m2cm}, {cm2m}},
			},
			"Energy loss due to ionizations": {
				DataItem: DataItem{
					saveFlag:   flag.Bool("li", false, "save ionization energy losses"),
					fileSuffix: "li",
				},
				columnNames: []string{"eV", "cm ^ -1"},
				values: func(de *DataExtractor) (args []float64, values [][]float64, labels []string) {
					for label := range de.energyLossByProcess {
						labels = append(labels, label)
					}
					sort.Strings(labels)
					for x := 0; x < de.model.numCells; x++ {
						args = append(args, de.model.xStep*float64(x))
						var row []float64
						for _, label := range labels {
							row = append(row, de.energyLossByProcess[label][x])
						}
						values = append(values, row)
					}
					return args, values, labels
				},
				scalers: [][]func(float64) float64{{m2cm}},
			},
		},
	}
}

func (de *DataExtractor) extract(modelName string, model *Model, parameters *ModelParameters) {
	if de.flags.outputPath != "" && de.flags.outputPath[len(de.flags.outputPath)-1] != '/' {
		de.flags.outputPath += "/"
	}

	de.model = model

	de.model.numCells = model.numCells

	de.crossSections = model.crossSections
	de.gasDensity = model.gasDensity

	de.cathodeFlux = parameters.CathodeCurrentDensity / electronCharge // [m^-2 s^-1]
	de.electronDensity = make([]float64, model.numCells)
	de.rateIntegral = make([]float64, model.numCells)
	de.sourceTerm = make([]float64, model.numCells)
	de.meanEnergy = make([]float64, model.numCells)
	de.driftVelocity = make([]float64, model.numCells)
	de.flux = make([]float64, model.numCells)
	de.velocity = make([]float64, model.numCells)
	de.radialVelocity = make([]float64, model.numCells)

	de.angularDistribution = make([][]float64, model.numCells)
	de.energyDistribution = make([][]float64, model.numCells)
	for x := range de.angularDistribution {
		de.angularDistribution[x] = make([]float64, model.numCellsMu)
		de.energyDistribution[x] = make([]float64, model.numCellsE)
	}

	de.TownsendAlphaF = make([]float64, model.numCells)
	de.collisions = make(map[string][]float64)
	de.energyLossByProcess = model.energyLossByProcess
	de.probabilitiesByProcess = make(map[string][]float64)
	de.probabilitiesByProcessFromEnergy = make(map[string][]float64)

	de.psiFIncrement = de.cathodeFlux / (float64(parameters.NElectrons) * parameters.DeltaE * parameters.DeltaMu)

	var lookUpVelocity []float64 = make([]float64, model.numCellsE+1)

	var energyRoot2Velocity float64 = math.Sqrt(2. / me)
	for eIndex := range lookUpVelocity {
		lookUpVelocity[eIndex] = math.Sqrt(eV2J(parameters.DeltaE*float64(eIndex)+parameters.DeltaE/2)) * energyRoot2Velocity
	}

	for xIndex := 0; xIndex < model.numCells; xIndex++ {
		for eIndex := 0; eIndex < model.numCellsE; eIndex++ {
			actualEnergy := model.parameters.DeltaE*float64(eIndex) + model.parameters.DeltaE/2.
			fXE := 0.
			fXE_mu := 0.
			v_x_fXE := 0.
			v_fXE := 0.
			v_r_fXE := 0.
			for muIndex := 0; muIndex < model.numCellsMu; muIndex++ {
				actualMu := model.parameters.DeltaMu*float64(muIndex) - 1.

				f := 0.
				if math.Abs(actualMu) > 0.0001 {
					f = de.psiFIncrement / (lookUpVelocity[eIndex] * math.Abs(actualMu)) * float64(model.distribution[xIndex][eIndex][muIndex])
				}
				fXE_mu += f * actualMu
				fXE += f

				v_x_fXE += f * lookUpVelocity[eIndex] * actualMu
				v_fXE += f * lookUpVelocity[eIndex]
				v_r_fXE += f * lookUpVelocity[eIndex] * math.Sqrt(1-actualMu*actualMu)
			}
			fXE_mu *= parameters.DeltaMu
			fXE *= parameters.DeltaMu
			v_x_fXE *= parameters.DeltaMu
			v_fXE *= parameters.DeltaMu
			v_r_fXE *= parameters.DeltaMu

			de.electronDensity[xIndex] += fXE

			de.meanEnergy[xIndex] += fXE * actualEnergy
			de.driftVelocity[xIndex] += fXE * lookUpVelocity[eIndex]
			de.flux[xIndex] += v_x_fXE
			de.velocity[xIndex] += v_fXE
			de.radialVelocity[xIndex] += v_r_fXE
			de.rateIntegral[xIndex] += model.crossSections.TotalCrossSectionOfKindAt(lxgata.IONIZATION, actualEnergy) * lookUpVelocity[eIndex] * fXE
		}
		de.electronDensity[xIndex] *= parameters.DeltaE

		de.meanEnergy[xIndex] *= parameters.DeltaE
		de.meanEnergy[xIndex] /= de.electronDensity[xIndex]

		de.driftVelocity[xIndex] *= parameters.DeltaE
		de.driftVelocity[xIndex] /= de.electronDensity[xIndex]

		de.flux[xIndex] *= parameters.DeltaE
		de.velocity[xIndex] *= parameters.DeltaE
		de.velocity[xIndex] /= de.electronDensity[xIndex]
		de.radialVelocity[xIndex] *= parameters.DeltaE
		de.radialVelocity[xIndex] /= de.electronDensity[xIndex]

		de.rateIntegral[xIndex] *= parameters.DeltaE
		de.rateIntegral[xIndex] /= de.electronDensity[xIndex]
		de.sourceTerm[xIndex] = model.gasDensity * de.electronDensity[xIndex] * de.rateIntegral[xIndex] // dn/dt = N * n(x) * <sigma_i(x,e) * v(x,e)>
	}

	de.model.xStep = model.xStep

	for processType, probs := range model.probabilities {
		de.probabilitiesByProcess[processType] = make([]float64, de.model.numCells)
		for xIndex := 0; xIndex < model.numCells; xIndex++ {
			de.probabilitiesByProcess[processType][xIndex] = probs[xIndex].avg()
		}
	}
	// for _, collision := range *de.crossSections {
	// 	de.probabilitiesByProcessFromEnergy[string(collision.Type)] = make([]float64, de.model.numCells)
	// 	for xIndex := 0; xIndex < model.numCells; xIndex++ {
	// 		totalCS := de.crossSections.TotalCrossSectionAt(de.meanEnergy[xIndex])
	// 		de.probabilitiesByProcessFromEnergy[string(collision.Type)][xIndex] += collision.CrossSectionAt(de.meanEnergy[xIndex]) / totalCS
	// 	}
	// }

	// initElectronDensity := model.cathodeFlux / math.Sqrt(2.*eV2J(4.5)/me)
	collCounters := make(map[string]float64)

	for key, val := range model.collisionAtCell {
		de.collisions[key] = make([]float64, de.model.numCells)
		for xIndex := 0; xIndex < model.numCells; xIndex++ {
			de.collisions[key][xIndex] = float64(val[xIndex]) / (float64(model.parameters.NElectrons) * model.xStep)
			collCounters[key] += float64(model.collisionAtCell[key][xIndex])
		}
	}
	for key := range collCounters {
		collCounters[key] /= float64(model.parameters.NElectrons)
	}
	fmt.Printf("Avg collisions per electron at distance %f: %v\n", model.parameters.CathodeFallLength, collCounters)

	for xIndex := 1; xIndex+2 < model.numCells; xIndex++ {
		de.TownsendAlphaF[xIndex] = (de.flux[xIndex+1] + de.flux[xIndex+2] - de.flux[xIndex] - de.flux[xIndex-1]) / (2. * de.flux[xIndex] * de.model.xStep)
	}
}

func openFile(makeDir bool, outputPath string, fileSuffix, modelName string) (*os.File, error) {
	if makeDir && fileSuffix != "" && fileSuffix != "." {
		os.MkdirAll(outputPath+fileSuffix, 0750)
		return os.Create(outputPath + fileSuffix + "/" + modelName + ".txt")
	} else {
		return os.Create(outputPath + modelName + "_" + fileSuffix + ".txt")
	}
}

func (de *DataExtractor) save(modelName string, parameters *ModelParameters) {
	for name, output := range de.flags.sequentials {
		if *output.saveFlag || *de.flags.all {
			var file *os.File
			file, err := openFile(parameters.MakeDir, de.flags.outputPath, output.fileSuffix, modelName)
			if err != nil {
				println("unable to save "+name+": ", err)
			} else {
				rows := [][]string{output.columnNames}
				xColumnValue, yColumnValues, yLabels := output.values(de)
				rows = append(rows, append([]string{""}, yLabels...))
				for x := range xColumnValue {
					var yColumnValuesStr []string
					if len(output.scalers) > 0 && output.scalers[0] != nil {
						for _, scaler := range output.scalers[0] {
							xColumnValue[x] = scaler(xColumnValue[x])
						}
					}
					for i := range yColumnValues[x] {
						if len(output.scalers) > 1 && output.scalers[1] != nil {
							for _, scaler := range output.scalers[1] {
								yColumnValues[x][i] = scaler(yColumnValues[x][i])
							}
						}
						yColumnValuesStr = append(yColumnValuesStr, strconv.FormatFloat(yColumnValues[x][i], 'f', -1, 64))
					}

					rows = append(rows, append([]string{strconv.FormatFloat(xColumnValue[x], 'f', -1, 64)}, yColumnValuesStr...))
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
