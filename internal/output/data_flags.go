package output

import (
	"encoding/csv"
	"encoding/json"
	"flag"
	"fmt"
	"log"
	"math"
	"os"
	"strconv"

	"github.com/wildstyl3r/lxgata"
	"github.com/wildstyl3r/stmc/internal/config"
	"github.com/wildstyl3r/stmc/internal/constants"
	"github.com/wildstyl3r/stmc/internal/extensions"
	"github.com/wildstyl3r/stmc/internal/model"
	"github.com/wildstyl3r/stmc/internal/utils"
)

type DataItem struct {
	SaveFlag *bool
	shortID  string
}

type Row struct {
	index, value, margin, upperMargin float64
}

type DataFile struct {
	metricsName          string
	indexName, valueName string
	data                 []Row
	confIntervals        bool
	twoSidedMargin       bool
}

type sequentialDataItem struct {
	DataItem
	values func(*model.Model) []DataFile
	xUnit  []config.UnitElement
	yUnit  []config.UnitElement
}

type customDataItem struct {
	DataItem
	prepare func(*model.Model) ([]byte, utils.ExportFileType)
}

type DataFlags struct {
	all         *bool
	Sequentials map[string]sequentialDataItem
	Custom      map[string]customDataItem
}

const (
	Potential                                   = "Potential"
	ElectricField                               = "ElectricField"
	ElectricFieldFromPotential                  = "ElectricFieldFromPotential"
	LfromV                                      = "PositionFromPotential"
	NormalizedCollisionCounters                 = "NormalizedCollisionCounters"
	NormalizedCollisionCountersFromDistribution = "NormalizedCollisionCountersD"
	NormalizedOutDischargeCollisionCounters     = "NormalizedOutDischargeCollisionCounters"
	NormalizedWallLoss                          = "NormalizedWallLoss"
	CollisionCounters                           = "CollisionCounters"
	DetailedCollisionCounters                   = "DetailedCollisionCounters"
	PlasmaDensity                               = "PlasmaDensity"
	MeanEnergy                                  = "MeanEnergy"
	NormalizedIonizationD                       = "NormalizedIonizationD"
	MeanVelocityX                               = "MeanVelocityX"
	MeanVelocityR                               = "MeanVelocityR"
	MeanVelocityAbs                             = "MeanVelocity"
	RawDistribution                             = "RawDistribution"
)

func NewDataFlags() DataFlags {
	return DataFlags{
		all: flag.Bool("all", false, "save every available metric"),
		Sequentials: map[string]sequentialDataItem{
			Potential: {
				DataItem: DataItem{
					SaveFlag: flag.Bool("p", false, "save potential"),
					shortID:  "V",
				},
				values: func(model *model.Model) []DataFile {
					dataFile := DataFile{
						metricsName: Potential,
						indexName:   fmt.Sprintf("x (%s)", model.Parameters.OutputUnit(config.Length)),
						valueName:   "g (V)",
						data:        make([]Row, model.NumCells),
					}
					for x := range model.NumCells {
						dataFile.data[x].index = model.XStep * float64(x)
						dataFile.data[x].value = model.VfromL(model.XStep * float64(x))
					}
					return []DataFile{dataFile}
				},
				xUnit: []config.UnitElement{{Class: config.Length, Power: 1}},
				yUnit: []config.UnitElement{},
			},
			ElectricField: {
				DataItem: DataItem{
					SaveFlag: flag.Bool("ef", false, "save Electric field"),
					shortID:  "Efield",
				},
				values: func(model *model.Model) []DataFile {
					dataFile := DataFile{
						metricsName: ElectricField,
						indexName:   fmt.Sprintf("x (%s)", model.Parameters.OutputUnit(config.Length)),
						valueName:   fmt.Sprintf("E (V/%s)", model.Parameters.OutputUnit(config.Length)),
						data:        make([]Row, model.NumCells),
					}
					for x := range model.NumCells {
						dataFile.data[x].index = model.XStep * float64(x)
						dataFile.data[x].value = model.EFieldFromL(model.XStep * float64(x))
					}
					return []DataFile{dataFile}
				},
				xUnit: []config.UnitElement{{Class: config.Length, Power: 1}},
				yUnit: []config.UnitElement{{Class: config.Length, Power: -1}},
			},
			ElectricFieldFromPotential: {
				DataItem: DataItem{
					SaveFlag: flag.Bool("efp", false, "save Electric field at potential"),
					shortID:  "efp",
				},
				values: func(model *model.Model) []DataFile {
					dataFile := DataFile{
						metricsName: ElectricFieldFromPotential,
						indexName:   "g(V)",
						valueName:   fmt.Sprintf("E (V/%s)", model.Parameters.OutputUnit(config.Length)),
					}
					for i := range 100 {
						v := -((model.Vc)*float64(100-i)/100. + model.Va)
						dataFile.data = append(dataFile.data, Row{index: v, value: model.EFieldFromPotential(v)})
					}
					for i := range 100 {
						v := -model.Va * float64(100-i) / 100.
						dataFile.data = append(dataFile.data, Row{index: v, value: model.EFieldFromPotential(v)})
					}
					return []DataFile{dataFile}
				},
				xUnit: []config.UnitElement{{Class: config.Length, Power: 1}},
				yUnit: []config.UnitElement{{Class: config.Length, Power: -1}},
			},
			LfromV: {
				DataItem: DataItem{
					SaveFlag: flag.Bool("lv", false, "save x from v"),
					shortID:  "lv",
				},
				values: func(model *model.Model) []DataFile {
					dataFile := DataFile{
						metricsName: LfromV,
						indexName:   "g(V)",
						valueName:   fmt.Sprintf("x (%s)", model.Parameters.OutputUnit(config.Length)),
					}
					for i := range 100 {
						v := -((model.Vc)*float64(100-i)/100. + model.Va)
						dataFile.data = append(dataFile.data, Row{index: v, value: model.LfromV(v)})
					}
					for i := range 100 {
						v := -model.Va * float64(100-i) / 100.
						dataFile.data = append(dataFile.data, Row{index: v, value: model.LfromV(v)})
					}
					return []DataFile{dataFile}
				},
				xUnit: []config.UnitElement{},
				yUnit: []config.UnitElement{{Class: config.Length, Power: 1}},
			},
			NormalizedCollisionCounters: {
				DataItem: DataItem{
					SaveFlag: flag.Bool("ncc", true, "save single-electron collision counters divided by pressure"),
					shortID:  "ncc",
				},
				values: func(model *model.Model) (files []DataFile) {
					collisions := model.GetMetrics(extensions.SingleElectronCollisionRateKey).(map[lxgata.CollisionType]utils.GriddedInterval)
					collisionsMargin := model.GetMetrics(extensions.SingleElectronCollisionRateMarginKey).(map[lxgata.CollisionType]utils.GriddedInterval)
					for label := range collisions {
						dataFile := DataFile{
							metricsName:   NormalizedCollisionCounters + "/" + string(label),
							indexName:     fmt.Sprintf("x (%s)", model.Parameters.OutputUnit(config.Length)),
							valueName:     fmt.Sprintf("N_i(%s^{-1} %s^{-1})", model.Parameters.OutputUnit(config.Length), model.Parameters.OutputUnit(config.Pressure)),
							data:          make([]Row, model.NumCells),
							confIntervals: true,
						}
						for x := range model.NumCells {
							dataFile.data[x].index = model.XStep * (float64(x) + 0.5)
							dataFile.data[x].value = collisions[lxgata.CollisionType(label)].Values[x] / model.Parameters.Pressure
							dataFile.data[x].margin = collisionsMargin[lxgata.CollisionType(label)].Values[x] / model.Parameters.Pressure
						}
						files = append(files, dataFile)
					}
					return files
				},
				xUnit: []config.UnitElement{{Class: config.Length, Power: 1}},
				yUnit: []config.UnitElement{{Class: config.Length, Power: -1}, {Class: config.Pressure, Power: -1}},
			},
			NormalizedCollisionCountersFromDistribution: {
				DataItem: DataItem{
					SaveFlag: flag.Bool("nccd", true, "save single-electron collision counters divided by pressure calculated from the EEDF"),
					shortID:  "nccd",
				},
				values: func(model *model.Model) (files []DataFile) {
					collisions := model.GetMetrics(extensions.SingleElectronCollisionRateKey).(map[lxgata.CollisionType]utils.GriddedInterval)
					collisionsMargin := model.GetMetrics(extensions.SingleElectronCollisionRateMarginKey).(map[lxgata.CollisionType]utils.GriddedInterval)
					for label := range collisions {
						dataFile := DataFile{
							metricsName:   NormalizedCollisionCounters + "/" + string(label),
							indexName:     fmt.Sprintf("x (%s)", model.Parameters.OutputUnit(config.Length)),
							valueName:     fmt.Sprintf("N_i(%s^{-1} %s^{-1})", model.Parameters.OutputUnit(config.Length), model.Parameters.OutputUnit(config.Pressure)),
							data:          make([]Row, model.NumCells),
							confIntervals: true,
						}
						for x := range model.NumCells {
							dataFile.data[x].index = model.XStep * (float64(x) + 0.5)
							dataFile.data[x].value = collisions[lxgata.CollisionType(label)].Values[x] / model.Parameters.Pressure
							dataFile.data[x].margin = collisionsMargin[lxgata.CollisionType(label)].Values[x] / model.Parameters.Pressure
						}
						files = append(files, dataFile)
					}
					return files
				},
				xUnit: []config.UnitElement{{Class: config.Length, Power: 1}},
				yUnit: []config.UnitElement{{Class: config.Length, Power: -1}, {Class: config.Pressure, Power: -1}},
			},
			CollisionCounters: {
				DataItem: DataItem{
					SaveFlag: flag.Bool("cc", true, "save collision counters"),
					shortID:  "cc",
				},
				values: func(model *model.Model) (files []DataFile) {
					collisions := model.GetMetrics(extensions.SingleElectronCollisionRateKey).(map[lxgata.CollisionType]utils.GriddedInterval)
					collisionsMargin := model.GetMetrics(extensions.SingleElectronCollisionRateMarginKey).(map[lxgata.CollisionType]utils.GriddedInterval)
					fluxAtCathode := model.Parameters.CathodeCurrentDensity * model.GetMetrics(extensions.CathodeElectronCurrentFractionKey).(float64) / constants.ElementaryCharge
					for label := range collisions {
						dataFile := DataFile{
							metricsName:   CollisionCounters + "/" + string(label),
							indexName:     fmt.Sprintf("x (%s)", model.Parameters.OutputUnit(config.Length)),
							valueName:     fmt.Sprintf("S_i(%s^{-3} %s^{-1})", model.Parameters.OutputUnit(config.Length), model.Parameters.OutputUnit(config.Time)),
							data:          make([]Row, model.NumCells),
							confIntervals: true,
						}
						for x := range model.NumCells {
							dataFile.data[x].index = model.XStep * (float64(x) + 0.5)
							dataFile.data[x].value = collisions[lxgata.CollisionType(label)].Values[x] * fluxAtCathode
							dataFile.data[x].margin = collisionsMargin[lxgata.CollisionType(label)].Values[x] * fluxAtCathode
						}
						files = append(files, dataFile)
					}
					return files
				},
				xUnit: []config.UnitElement{{Class: config.Length, Power: 1}},
				yUnit: []config.UnitElement{{Class: config.Length, Power: -3}, {Class: config.Time, Power: -1}},
			},
			DetailedCollisionCounters: {
				DataItem: DataItem{
					SaveFlag: flag.Bool("dcc", true, "save detailed collision counters"),
					shortID:  "dcc",
				},
				values: func(model *model.Model) (files []DataFile) {
					collisions := model.GetMetrics(extensions.SingleElectronDetailedCollisionRateKey).(map[string]utils.GriddedInterval)
					// collisionsLowerMargin := model.GetMetrics(extensions.SingleElectronDetailedCollisionRateLowerMarginKey).(map[string]utils.GriddedInterval)
					// collisionsUpperMargin := model.GetMetrics(extensions.SingleElectronDetailedCollisionRateUpperMarginKey).(map[string]utils.GriddedInterval)
					fluxAtCathdode := model.Parameters.CathodeCurrentDensity * model.GetMetrics(extensions.CathodeElectronCurrentFractionKey).(float64) / constants.ElementaryCharge
					for label := range collisions {
						dataFile := DataFile{
							metricsName:   DetailedCollisionCounters + "/" + string(label),
							indexName:     fmt.Sprintf("x (%s)", model.Parameters.OutputUnit(config.Length)),
							valueName:     fmt.Sprintf("S_i(%s^{-3} %s^{-1})", model.Parameters.OutputUnit(config.Length), model.Parameters.OutputUnit(config.Time)),
							data:          make([]Row, model.NumCells),
							confIntervals: false,
						}
						for x := range model.NumCells {
							dataFile.data[x].index = model.XStep * (float64(x) + 0.5)
							dataFile.data[x].value = collisions[label].Values[x] * fluxAtCathdode
							// dataFile.data[x].margin = collisionsLowerMargin[label].Values[x] * fluxAtCathdode
							// dataFile.data[x].upperMargin = collisionsUpperMargin[label].Values[x] * fluxAtCathdode
						}
						files = append(files, dataFile)
					}
					return files
				},
				xUnit: []config.UnitElement{{Class: config.Length, Power: 1}},
				yUnit: []config.UnitElement{{Class: config.Length, Power: -3}, {Class: config.Time, Power: -1}},
			},
			NormalizedWallLoss: {
				DataItem: DataItem{
					SaveFlag: flag.Bool("nwl", true, "save normalized wall loss"),
					shortID:  "nwl",
				},
				values: func(model *model.Model) []DataFile {
					dataFile := DataFile{
						metricsName:   NormalizedWallLoss,
						indexName:     fmt.Sprintf("x (%s)", model.Parameters.OutputUnit(config.Length)),
						valueName:     fmt.Sprintf("N_i(%s^{-1} %s^{-1})", model.Parameters.OutputUnit(config.Length), model.Parameters.OutputUnit(config.Pressure)),
						data:          make([]Row, model.NumCells),
						confIntervals: true,
					}
					wallLosses := model.GetMetrics(extensions.NormalizedWallLossRateKey).([]float64)
					wallLossesMargin := model.GetMetrics(extensions.NormalizedWallLossRateMarginKey).([]float64)
					for x := range model.NumCells {
						dataFile.data[x].index = model.XStep * (float64(x) + 0.5)
						dataFile.data[x].value = wallLosses[x]
						dataFile.data[x].margin = wallLossesMargin[x]
					}
					return []DataFile{dataFile}
				},
				xUnit: []config.UnitElement{{Class: config.Length, Power: 1}},
				yUnit: []config.UnitElement{{Class: config.Length, Power: -1}, {Class: config.Pressure, Power: -1}},
			},
			NormalizedOutDischargeCollisionCounters: {
				DataItem: DataItem{
					SaveFlag: flag.Bool("nocc", true, "save single-electron out of discharge volume collision counters divided by pressure"),
					shortID:  "nocc",
				},
				values: func(model *model.Model) (files []DataFile) {
					collisions := model.GetMetrics(extensions.SingleElectronCollisionsOutsideKey).(map[lxgata.CollisionType]utils.GriddedInterval)
					for label := range collisions {
						dataFile := DataFile{
							metricsName:   NormalizedOutDischargeCollisionCounters + "/" + string(label),
							indexName:     fmt.Sprintf("x (%s)", model.Parameters.OutputUnit(config.Length)),
							valueName:     fmt.Sprintf("N_i(%s^{-1} %s^{-1})", model.Parameters.OutputUnit(config.Length), model.Parameters.OutputUnit(config.Pressure)),
							data:          make([]Row, model.NumCells),
							confIntervals: true,
						}
						for x := range model.NumCells {
							dataFile.data[x].index = model.XStep * (float64(x) + 0.5)
							dataFile.data[x].value = collisions[lxgata.CollisionType(label)].Values[x] / model.Parameters.Pressure
						}
						files = append(files, dataFile)
					}
					return files
				},
				xUnit: []config.UnitElement{{Class: config.Length, Power: 1}},
				yUnit: []config.UnitElement{{Class: config.Length, Power: -1}, {Class: config.Pressure, Power: -1}},
			},
			// "Energy loss due to ionizations": {
			// 	DataItem: DataItem{
			// 		SaveFlag:   flag.Bool("li", false, "save ionization energy losses"),
			// 		fileSuffix: "li",
			// 	},
			// 	values: func(model *model.Model) []DataFile {
			// 		for x := range model.NumCells {
			// 			args = append(args, model.XStep*(float64(x)+0.5))
			// 			row := []float64{model.EnergyLossByProcess[lxgata.IONIZATION][x]}
			// 			values = append(values, row)
			// 		}
			// 		return args, values, []string{string(lxgata.IONIZATION)}, []string{
			// 			fmt.Sprintf("x (%s)", model.Parameters.OutputUnit(config.Length)),
			// 			fmt.Sprintf("(%s^{-1}%s^{-1})", model.Parameters.OutputUnit(config.Energy), model.Parameters.OutputUnit(config.Length))}
			// 	},
			// 	xUnit: []config.UnitElement{{Class: config.Length, Power: 1}},
			// 	yUnit: []config.UnitElement{{Class: config.Length, Power: -1}},
			// },
			// "Out of energy for ionizations": {
			// 	DataItem: DataItem{
			// 		SaveFlag:   flag.Bool("oe", false, "save ooe exit events count"),
			// 		fileSuffix: "oe",
			// 	},
			// 	values: func(model *model.Model) []DataFile {
			// 		collisions := model.GetMetrics(extensions.NormalizedCollisionRateKey).(map[lxgata.CollisionType][]float64)
			// 		for label := range collisions {
			// 			labels = append(labels, string(label))
			// 			labels = append(labels, string(label)+"_conf_interval")
			// 		}
			// 		sort.Strings(labels)
			// 		for x := range model.NumCells {
			// 			args = append(args, model.XStep*(float64(x)+0.5))
			// 			var row []float64
			// 			row = append(row, float64((model.OutOfEnergyAtCell[x]))/float64(model.Parameters.NElectrons))
			// 			values = append(values, row)
			// 		}
			// 		return args, values, nil, []string{
			// 			fmt.Sprintf("x (%s)", model.Parameters.OutputUnit(config.Length)),
			// 			fmt.Sprintf("(%s^{-1}%s^{-1})", model.Parameters.OutputUnit(config.Energy), model.Parameters.OutputUnit(config.Length))}
			// 	},
			// 	xUnit: []config.UnitElement{{Class: config.Length, Power: 1}},
			// 	yUnit: []config.UnitElement{{Class: config.Length, Power: -1}},
			// },
			PlasmaDensity: {
				DataItem: DataItem{
					SaveFlag: flag.Bool("n", true, "save plasma density"),
					shortID:  "n",
				},
				values: func(model *model.Model) []DataFile {
					dataFile := DataFile{
						metricsName: PlasmaDensity,
						indexName:   fmt.Sprintf("x (%s)", model.Parameters.OutputUnit(config.Length)),
						valueName:   fmt.Sprintf("(%s^{-3})", model.Parameters.OutputUnit(config.Length)),
						data:        make([]Row, model.NumCells),
					}
					densityPerElectronFlux := model.GetMetrics(extensions.DensityPerCathodeElectronFluxKey).(utils.GriddedInterval)
					fluxAtCathode := model.Parameters.CathodeCurrentDensity * model.GetMetrics(extensions.CathodeElectronCurrentFractionKey).(float64) / constants.ElementaryCharge

					for x := range model.NumCells {
						dataFile.data[x].index = model.XStep * (float64(x) + 0.5)
						dataFile.data[x].value = densityPerElectronFlux.Values[x] * fluxAtCathode
					}
					return []DataFile{dataFile}
				},
				xUnit: []config.UnitElement{{Class: config.Length, Power: 1}},
				yUnit: []config.UnitElement{{Class: config.Length, Power: -3}},
			},
			MeanEnergy: {
				DataItem: DataItem{
					SaveFlag: flag.Bool("e", false, "save mean energy"),
					shortID:  "e",
				},
				values: func(model *model.Model) []DataFile {
					dataFile := DataFile{
						metricsName: MeanEnergy,
						indexName:   fmt.Sprintf("x (%s)", model.Parameters.OutputUnit(config.Length)),
						valueName:   fmt.Sprintf("\\varepsilon(%s)", model.Parameters.OutputUnit(config.Energy)),
						data:        make([]Row, model.NumCells),
					}

					electronDensity := make([]float64, model.NumCells)
					meanEnergy := make([]float64, model.NumCells)

					var lookUpVelocity []float64 = make([]float64, model.NumCellsE+1)

					var energyRoot2Velocity float64 = math.Sqrt(2. / constants.ElectornMass)
					for eIndex := range lookUpVelocity {
						lookUpVelocity[eIndex] = math.Sqrt(utils.EV2J(model.Parameters.EnergyDiscretizationStep*(float64(eIndex)+0.5))) * energyRoot2Velocity
					}

					for xIndex := 0; xIndex < model.NumCells; xIndex++ {
						for eIndex := 0; eIndex < model.NumCellsE; eIndex++ {
							currentEnergy := model.Parameters.EnergyDiscretizationStep * (float64(eIndex) + 0.5)
							fXE := 0.
							fXECompensation := 0.
							for muIndex := 0; muIndex < model.NumCellsMu; muIndex++ {
								currentMu := max(-1, min(model.Parameters.MuDiscretizationStep*(float64(muIndex)+0.5)-1., 1))

								f := 0.
								if math.Abs(currentMu) > 0.0001 {
									f = 1 / (lookUpVelocity[eIndex] * math.Abs(currentMu)) * float64(model.DistributionXEMu[xIndex][eIndex][muIndex])
								}
								y := f - fXECompensation
								temp := fXE + y
								fXECompensation = (temp - fXE) - y
								fXE = temp
							}

							electronDensity[xIndex] += fXE

							meanEnergy[xIndex] += fXE * currentEnergy
						}
						meanEnergy[xIndex] /= electronDensity[xIndex]
						dataFile.data[xIndex].index = model.XStep * (float64(xIndex) + 0.5)
						dataFile.data[xIndex].value = meanEnergy[xIndex]
					}
					return []DataFile{dataFile}
				},
				xUnit: []config.UnitElement{{Class: config.Length, Power: 1}},
				yUnit: []config.UnitElement{{Class: config.Energy, Power: 1}},
			},
			NormalizedIonizationD: {
				DataItem: DataItem{
					SaveFlag: flag.Bool("nicd", true, "save normalized source term calculated from EEDF"),
					shortID:  "nicd",
				},
				values: func(model *model.Model) []DataFile {
					dataFile := DataFile{
						metricsName: NormalizedIonizationD,
						indexName:   fmt.Sprintf("x (%s)", model.Parameters.OutputUnit(config.Length)),
						valueName:   fmt.Sprintf("S (%s^-1 %s^-1)", model.Parameters.OutputUnit(config.Length), model.Parameters.OutputUnit(config.Pressure)),
						data:        make([]Row, model.NumCells),
					}

					// ionSource := make([]float64, model.NumCells)

					var lookUpVelocity []float64 = make([]float64, model.NumCellsE+1)

					// var energyRoot2Velocity float64 = math.Sqrt(2. * constants.ElectronChargeToMassRatio)
					for eIndex := range lookUpVelocity {
						energy := model.Parameters.EnergyDiscretizationStep * (float64(eIndex) + 0.5)
						lookUpVelocity[eIndex] = utils.EV2electronVelocity(energy)
					}

					var ionizationCS *lxgata.Collision
					for p := range model.Parameters.CrossSectionsData().Processes {
						if model.Parameters.CrossSectionsData().Processes[p].Type == lxgata.IONIZATION {
							ionizationCS = &model.Parameters.CrossSectionsData().Processes[p]
							break
						}
					}
					if ionizationCS == nil {
						panic("no ionization cs found")
					}

					for x := range model.DistributionXEMu {
						norm, normCompensation := 0., 0.
						ionSourceX, ionSourceXCompensation := 0., 0.

						for e := range model.DistributionXEMu[x] {
							currentEnergy := model.Parameters.EnergyDiscretizationStep * (float64(e) + 0.5)
							fXE, fXECompensation := 0., 0.
							for mu := range model.DistributionXEMu[x][e] {
								utils.KahanSum(&fXE, &model.DistributionXEMu[x][e][mu], &fXECompensation)
							}
							fXE *= model.Parameters.MuDiscretizationStep
							normIntegrable := math.Sqrt(currentEnergy) * fXE * model.Parameters.EnergyDiscretizationStep
							utils.KahanSum(&norm, &normIntegrable, &normCompensation)
							sourceIntegrable := utils.EV2electronVelocity(currentEnergy) * ionizationCS.CrossSectionAt(currentEnergy) * normIntegrable
							utils.KahanSum(&ionSourceX, &sourceIntegrable, &ionSourceXCompensation)

							//norm /= float64(model.TotalElectronsEmittedOnCathode)
						}
						fmt.Printf("norm: %.9f\n", norm)
						dataFile.data[x].index = model.XStep * (float64(x) + 0.5)
						dataFile.data[x].value = ionSourceX / norm / model.Parameters.Pressure
					}

					// for xIndex := 0; xIndex < model.NumCells; xIndex++ {
					// 	ionSourceXCompensation := 0.
					// 	for eIndex := int(ionizationCS.Threshold / model.Parameters.EnergyDiscretizationStep); eIndex < model.NumCellsE; eIndex++ {
					// 		currentEnergy := model.Parameters.EnergyDiscretizationStep * (float64(eIndex) + 0.5)
					// 		fXE := 0.
					// 		fXECompensation := 0.
					// 		for muIndex := 0; muIndex < model.NumCellsMu; muIndex++ {
					// 			currentMu := max(-1, min(model.Parameters.MuDiscretizationStep*(float64(muIndex)+0.5)-1., 1))

					// 			f := 0.
					// 			if math.Abs(currentMu) > 0.00001 {
					// 				f = 1. / (lookUpVelocity[eIndex] * math.Abs(currentMu) * float64(model.TotalElectronsEmittedOnCathode) * model.Parameters.MuDiscretizationStep * model.Parameters.EnergyDiscretizationStep) * float64(model.DistributionXEMu[xIndex][eIndex][muIndex])
					// 			}
					// 			y := f - fXECompensation
					// 			temp := fXE + y
					// 			fXECompensation = (temp - fXE) - y
					// 			fXE = temp
					// 		}
					// 		y := fXE*currentEnergy*ionizationCS.CrossSectionAt(currentEnergy) - ionSourceXCompensation
					// 		temp := ionSource[xIndex] + y
					// 		ionSourceXCompensation = (temp - ionSource[xIndex]) - y
					// 		ionSource[xIndex] = temp
					// 	}
					// 	dataFile.data[xIndex].index = model.XStep * (float64(xIndex) + 0.5)
					// 	dataFile.data[xIndex].value = ionSource[xIndex] * energyRoot2Velocity / model.Parameters.Pressure
					// }
					return []DataFile{dataFile}
				},
				xUnit: []config.UnitElement{{Class: config.Length, Power: 1}},
				yUnit: []config.UnitElement{{Class: config.Length, Power: -1}, {Class: config.Pressure, Power: -1}},
			},
			MeanVelocityX: {
				DataItem: DataItem{
					SaveFlag: flag.Bool("vx", false, "save drift velocity"),
					shortID:  "vx",
				},
				values: func(model *model.Model) []DataFile {
					dataFile := DataFile{
						metricsName: MeanVelocityX,
						indexName:   fmt.Sprintf("x (%s)", model.Parameters.OutputUnit(config.Length)),
						valueName:   fmt.Sprintf("v_x(%s%s^{-1})", model.Parameters.OutputUnit(config.Length), model.Parameters.OutputUnit(config.Time)),
						data:        make([]Row, model.NumCells),
					}

					electronDensity := make([]float64, model.NumCells)
					flux := make([]float64, model.NumCells)

					var lookUpVelocity []float64 = make([]float64, model.NumCellsE+1)

					var energyRoot2Velocity float64 = math.Sqrt(2. / constants.ElectornMass)
					for eIndex := range lookUpVelocity {
						lookUpVelocity[eIndex] = math.Sqrt(utils.EV2J(model.Parameters.EnergyDiscretizationStep*(float64(eIndex)+0.5))) * energyRoot2Velocity
					}

					for xIndex := range model.NumCells {
						for eIndex := range model.NumCellsE {
							fXE := 0.
							v_x_fXE := 0.
							for muIndex := range model.NumCellsMu {
								currentMu := max(-1, min(model.Parameters.MuDiscretizationStep*(float64(muIndex)+0.5)-1., 1))

								f := 0.
								if math.Abs(currentMu) > 0.0001 {
									f = 1 / (lookUpVelocity[eIndex] * math.Abs(currentMu)) * float64(model.DistributionXEMu[xIndex][eIndex][muIndex])
								}
								fXE += f

								v_x_fXE += f * lookUpVelocity[eIndex] * currentMu
							}

							electronDensity[xIndex] += fXE
							flux[xIndex] += v_x_fXE
						}
						dataFile.data[xIndex].index = model.XStep * (float64(xIndex) + 0.5)
						dataFile.data[xIndex].value = flux[xIndex] / electronDensity[xIndex]
					}
					return []DataFile{dataFile}
				},
				xUnit: []config.UnitElement{{Class: config.Length, Power: 1}},
				yUnit: []config.UnitElement{{Class: config.Length, Power: 1}, {Class: config.Time, Power: -1}},
			},
		},
		Custom: map[string]customDataItem{
			RawDistribution: {
				DataItem: DataItem{
					SaveFlag: flag.Bool("raw", false, "save raw distribution"),
					shortID:  "raw",
				},
				prepare: func(m *model.Model) ([]byte, utils.ExportFileType) {
					jsonData, err := json.Marshal(m.DistributionXEMu)
					if err != nil {
						print("unable to jsonify raw distribution", err)
						return nil, utils.ExportFileType("")
					} else {
						return jsonData, utils.TypeJson
					}
				},
			},
		},
	}
}

func Save(modelName string, model *model.Model, df DataFlags, outputPath string) {

	for _, output := range df.Sequentials {
		if *output.SaveFlag || *df.all {
			drafts := output.values(model)
			for f := range drafts {
				path := outputPath + modelName
				file, err := utils.OpenFile(path+"/"+drafts[f].metricsName, utils.TypeCSV)
				if err != nil {
					fmt.Printf("unable to save %s: %+v\n", path+"/"+drafts[f].metricsName, err)
				} else {
					rows := [][]string{
						{drafts[f].indexName, drafts[f].valueName},
					}
					if drafts[f].confIntervals {
						if drafts[f].twoSidedMargin {
							rows[0] = append(rows[0], "marginL", "marginU")
							for i := range drafts[f].data {
								rows = append(rows, []string{
									strconv.FormatFloat(config.FieldToSIeV(drafts[f].data[i].index, output.xUnit, model.Parameters.OutputUnits(), false), 'f', -1, 64),
									strconv.FormatFloat(config.FieldToSIeV(drafts[f].data[i].value, output.yUnit, model.Parameters.OutputUnits(), false), 'f', -1, 64),
									strconv.FormatFloat(config.FieldToSIeV(drafts[f].data[i].margin, output.yUnit, model.Parameters.OutputUnits(), false), 'f', -1, 64),
									strconv.FormatFloat(config.FieldToSIeV(drafts[f].data[i].upperMargin, output.yUnit, model.Parameters.OutputUnits(), false), 'f', -1, 64),
								})

							}
						} else {
							rows[0] = append(rows[0], "margin")
							for i := range drafts[f].data {
								rows = append(rows, []string{
									strconv.FormatFloat(config.FieldToSIeV(drafts[f].data[i].index, output.xUnit, model.Parameters.OutputUnits(), false), 'f', -1, 64),
									strconv.FormatFloat(config.FieldToSIeV(drafts[f].data[i].value, output.yUnit, model.Parameters.OutputUnits(), false), 'f', -1, 64),
									strconv.FormatFloat(config.FieldToSIeV(drafts[f].data[i].margin, output.yUnit, model.Parameters.OutputUnits(), false), 'f', -1, 64),
								})

							}
						}

					} else {
						for i := range drafts[f].data {
							rows = append(rows, []string{
								strconv.FormatFloat(config.FieldToSIeV(drafts[f].data[i].index, output.xUnit, model.Parameters.OutputUnits(), false), 'f', -1, 64),
								strconv.FormatFloat(config.FieldToSIeV(drafts[f].data[i].value, output.yUnit, model.Parameters.OutputUnits(), false), 'f', -1, 64),
							})
						}
					}
					w := csv.NewWriter(file)
					w.WriteAll(rows)
					if model.Parameters.Verbose() {
						print(output.shortID, " ")
					}
					if err := w.Error(); err != nil {
						log.Fatalln("error writing csv:", err)
					}
				}
			}
		}
	}
	println("")
	for name, output := range df.Custom {
		if *output.SaveFlag || *df.all {
			data, fileExt := output.prepare(model)
			if data != nil {
				var file *os.File
				file, err := utils.OpenFile(outputPath+"/"+modelName+"/"+output.shortID, fileExt)
				if err != nil {
					println("unable to save "+name+": ", err)
				} else {
					file.Write(data)
					file.Sync()
					file.Close()
					if model.Parameters.Verbose() {
						println(name + " saved")
					}
				}
			}

		}
	}
}
