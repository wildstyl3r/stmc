package output

import (
	"encoding/csv"
	"flag"
	"log"
	"math"
	"os"
	"sort"
	"strconv"

	"github.com/wildstyl3r/lxgata"
	"github.com/wildstyl3r/stmc/internal/config"
	"github.com/wildstyl3r/stmc/internal/datahub"
	"github.com/wildstyl3r/stmc/internal/extensions"
	"github.com/wildstyl3r/stmc/internal/model"
	"github.com/wildstyl3r/stmc/internal/utils"
)

type DataItem struct {
	SaveFlag   *bool
	fileSuffix string
}

type sequentialDataItem struct {
	DataItem
	columnNames []string
	values      func(*model.Model) (args []float64, values [][]float64, labels []string)
	xUnit       []config.UnitElement
	yUnit       []config.UnitElement
}

type DataFlags struct {
	all         *bool
	Sequentials map[string]sequentialDataItem
}

const (
	Potential                  = "Potential"
	ElectricField              = "Electric field"
	ElectricFieldFromPotential = "Electric field from potential"
	LfromV                     = "Position from potential"
	CollisionCounters          = "Collision counters"
	NormalizedSourceTerm       = "Normalized source term"
	PlasmaDensity              = "Plasma density"
	MeanEnergy                 = "Mean energy"
	MeanVelocityX              = "Mean velocity along x"
)

func NewDataFlags() DataFlags {
	return DataFlags{
		all: flag.Bool("all", false, "save every available metric"),
		Sequentials: map[string]sequentialDataItem{
			Potential: {
				DataItem: DataItem{
					SaveFlag:   flag.Bool("p", false, "save potential"),
					fileSuffix: "V",
				},
				columnNames: []string{"x (cm)", "g (V)"},
				values: func(model *model.Model) (args []float64, values [][]float64, labels []string) {
					for x := range model.NumCells {
						args = append(args, model.XStep*float64(x))
						values = append(values, []float64{model.VfromL(model.XStep * float64(x))})
					}
					return args, values, nil
				},
				xUnit: []config.UnitElement{{Class: config.Length, Power: 1}},
				yUnit: []config.UnitElement{},
			},
			ElectricField: {
				DataItem: DataItem{
					SaveFlag:   flag.Bool("ef", false, "save Electric field"),
					fileSuffix: "Efield",
				},
				columnNames: []string{"x (cm)", "E (V/m)"},
				values: func(model *model.Model) (args []float64, values [][]float64, labels []string) {
					for x := range model.NumCells {
						args = append(args, model.XStep*float64(x))
						values = append(values, []float64{model.EFieldFromL(model.XStep * float64(x))})
					}
					return args, values, nil
				},
				xUnit: []config.UnitElement{{Class: config.Length, Power: 1}},
				yUnit: []config.UnitElement{{Class: config.Length, Power: -1}},
			},
			ElectricFieldFromPotential: {
				DataItem: DataItem{
					SaveFlag:   flag.Bool("efp", false, "save Electric field at potential"),
					fileSuffix: "efp",
				},
				columnNames: []string{"V", "E (V/m)"},
				values: func(model *model.Model) (args []float64, values [][]float64, labels []string) {
					for i := range 100 {
						v := -((model.Vc)*float64(100-i)/100. + model.Va)
						args = append(args, v)
						values = append(values, []float64{model.EFieldFromPotential(v)})
					}
					for i := range 100 {
						v := -model.Va * float64(100-i) / 100.
						args = append(args, v)
						values = append(values, []float64{model.EFieldFromPotential(v)})
					}
					return args, values, nil
				},
				xUnit: []config.UnitElement{{Class: config.Length, Power: 1}},
				yUnit: []config.UnitElement{{Class: config.Length, Power: -1}},
			},
			LfromV: {
				DataItem: DataItem{
					SaveFlag:   flag.Bool("lv", false, "save x from v"),
					fileSuffix: "lv",
				},
				columnNames: []string{"g (V)", "x (cm)"},
				values: func(model *model.Model) (args []float64, values [][]float64, labels []string) {
					for i := range 100 {
						v := -((model.Vc)*float64(100-i)/100. + model.Va)
						args = append(args, v)
						values = append(values, []float64{model.LfromV(v)})
					}
					for i := range 100 {
						v := -model.Va * float64(100-i) / 100.
						args = append(args, v)
						values = append(values, []float64{model.LfromV(v)})
					}
					return args, values, nil
				},
				xUnit: []config.UnitElement{},
				yUnit: []config.UnitElement{{Class: config.Length, Power: 1}},
			},
			CollisionCounters: {
				DataItem: DataItem{
					SaveFlag:   flag.Bool("cc", false, "save collision counters"),
					fileSuffix: "cc",
				},
				columnNames: []string{"x (cm)", "N_i(cm^{-1} Torr^{-1})", "Standard error"},
				values: func(model *model.Model) (args []float64, values [][]float64, labels []string) {
					collisions := datahub.Get(extensions.NormalizedCollisionRateKey, model).(map[lxgata.CollisionType][]float64)
					collisionsCI := datahub.Get(extensions.NormalizedCollisionRateCIKey, model).(map[lxgata.CollisionType][]float64)
					for label := range collisions {
						labels = append(labels, string(label))
						if model.Parameters.CalculateStdError {
							labels = append(labels, string(label)+"_conf_interval")
						}
					}
					sort.Strings(labels)
					for x := range model.NumCells {
						args = append(args, model.XStep*(float64(x)+0.5))
						var row []float64
						for _, label := range labels {
							row = append(row, collisions[lxgata.CollisionType(label)][x])
							if model.Parameters.CalculateStdError {
								row = append(row, collisionsCI[lxgata.CollisionType(label)][x])
							}
						}
						values = append(values, row)
					}
					return args, values, labels
				},
				xUnit: []config.UnitElement{{Class: config.Length, Power: 1}},
				yUnit: []config.UnitElement{{Class: config.Length, Power: -1}, {Class: config.Pressure, Power: -1}},
			},
			NormalizedSourceTerm: {
				DataItem: DataItem{
					SaveFlag:   flag.Bool("nst", true, "save normalized source term"),
					fileSuffix: "nst",
				},
				columnNames: []string{"x (cm)", "N_i(cm^{-1} Torr^{-1})", "Standard error"},
				values: func(model *model.Model) (args []float64, values [][]float64, labels []string) {
					collisions := datahub.Get(extensions.NormalizedCollisionRateKey, model).(map[lxgata.CollisionType][]float64)
					collisionsCI := datahub.Get(extensions.NormalizedCollisionRateCIKey, model).(map[lxgata.CollisionType][]float64)
					labels = append(labels, string(lxgata.IONIZATION))
					if model.Parameters.CalculateStdError {
						labels = append(labels, string(lxgata.IONIZATION)+"_conf_interval")
					}
					sort.Strings(labels)
					for x := range model.NumCells {
						args = append(args, model.XStep*(float64(x)+0.5))
						var row []float64
						row = append(row, collisions[lxgata.IONIZATION][x])
						if model.Parameters.CalculateStdError {
							row = append(row, collisionsCI[lxgata.IONIZATION][x])
						}
						values = append(values, row)
					}
					return args, values, labels
				},
				xUnit: []config.UnitElement{{Class: config.Length, Power: 1}},
				yUnit: []config.UnitElement{{Class: config.Length, Power: -1}, {Class: config.Pressure, Power: -1}},
			},
			"Energy loss due to ionizations": {
				DataItem: DataItem{
					SaveFlag:   flag.Bool("li", false, "save ionization energy losses"),
					fileSuffix: "li",
				},
				columnNames: []string{"eV", "cm ^ -1"},
				values: func(model *model.Model) (args []float64, values [][]float64, labels []string) {
					for x := range model.NumCells {
						args = append(args, model.XStep*(float64(x)+0.5))
						row := []float64{model.EnergyLossByProcess[lxgata.IONIZATION][x]}
						values = append(values, row)
					}
					return args, values, []string{string(lxgata.IONIZATION)}
				},
				xUnit: []config.UnitElement{{Class: config.Length, Power: 1}},
				yUnit: []config.UnitElement{{Class: config.Length, Power: -1}},
			},
			"Out of energy for ionizations": {
				DataItem: DataItem{
					SaveFlag:   flag.Bool("oe", false, "save ooe exit events count"),
					fileSuffix: "oe",
				},
				columnNames: []string{"x (cm)", "cm ^ -1"},
				values: func(model *model.Model) (args []float64, values [][]float64, labels []string) {
					collisions := datahub.Get(extensions.NormalizedCollisionRateKey, model).(map[lxgata.CollisionType][]float64)
					for label := range collisions {
						labels = append(labels, string(label))
						if model.Parameters.CalculateStdError {
							labels = append(labels, string(label)+"_conf_interval")
						}
					}
					sort.Strings(labels)
					for x := range model.NumCells {
						args = append(args, model.XStep*(float64(x)+0.5))
						var row []float64
						row = append(row, float64((model.OutOfEnergyAtCell[x]))/float64(model.Parameters.NElectrons))
						values = append(values, row)
					}
					return args, values, nil
				},
				xUnit: []config.UnitElement{{Class: config.Length, Power: 1}},
				yUnit: []config.UnitElement{{Class: config.Length, Power: -1}},
			},
			PlasmaDensity: {
				DataItem: DataItem{
					SaveFlag:   flag.Bool("n", true, "save plasma density"),
					fileSuffix: "n",
				},
				columnNames: []string{"x (cm)", "cm ^ -3"},
				values: func(model *model.Model) (args []float64, values [][]float64, labels []string) {
					density := datahub.Get("GlowDischargeDensity", model).([]float64)
					for x := range model.NumCells {
						args = append(args, model.XStep*(float64(x)+0.5))
						values = append(values, []float64{density[x]})
					}
					return args, values, []string{"Plasma density n(x)"}
				},
				xUnit: []config.UnitElement{{Class: config.Length, Power: 1}},
				yUnit: []config.UnitElement{{Class: config.Length, Power: -3}},
			},
			MeanEnergy: {
				DataItem: DataItem{
					SaveFlag:   flag.Bool("e", false, "save mean energy"),
					fileSuffix: "e",
				},
				columnNames: []string{"x (cm)", "eV"},
				values: func(m *model.Model) (args []float64, values [][]float64, labels []string) {
					for x := range m.NumCells {
						meanEnergy := 0.
						n := 0.
						for j := range m.NumMuCells {
							mu := (float64(j-m.NumMuCells/2) + 0.5) * m.MuStep
							for iterator := m.Distribution[x][j].Front(); iterator != nil; iterator = iterator.Next() {
								energy := iterator.Value.(float64)
								divBy := mu * utils.EV2electronVelocity((math.Floor(energy/m.EStep)+0.5)*m.EStep)
								meanEnergy += energy / divBy
								n += 1. / divBy
							}
						}
						args = append(args, m.XStep*(float64(x)+0.5))
						values = append(values, []float64{meanEnergy / n})
					}
					return args, values, []string{"varepsilon (x)"}
				},
				xUnit: []config.UnitElement{{Class: config.Length, Power: 1}},
				yUnit: []config.UnitElement{{Class: config.Energy, Power: 1}},
			},
			MeanVelocityX: {
				DataItem: DataItem{
					SaveFlag:   flag.Bool("vx", false, "save drift velocity"),
					fileSuffix: "vx",
				},
				columnNames: []string{"x (cm)", "cm s^-1"},
				values: func(m *model.Model) (args []float64, values [][]float64, labels []string) {
					for x := range m.NumCells {
						vx := 0.
						n := 0.
						for j := range m.NumMuCells {
							mu := (float64(j-m.NumMuCells/2) + 0.5) * m.MuStep
							for iterator := m.Distribution[x][j].Front(); iterator != nil; iterator = iterator.Next() {
								energy := iterator.Value.(float64)
								divBy := mu * utils.EV2electronVelocity((math.Floor(energy/m.EStep)+0.5)*m.EStep)
								vx += 1. / mu
								n += 1. / divBy
							}
						}
						args = append(args, m.XStep*(float64(x)+0.5))
						values = append(values, []float64{vx / n})
					}
					return args, values, []string{"Vx (x)"}
				},
				xUnit: []config.UnitElement{{Class: config.Length, Power: 1}},
				yUnit: []config.UnitElement{{Class: config.Length, Power: 1}, {Class: config.Time, Power: -1}},
			},
		},
	}
}

func Save(modelName string, model *model.Model, df DataFlags, outputPath string) {
	for name, output := range df.Sequentials {
		if *output.SaveFlag || *df.all {
			var file *os.File
			file, err := utils.OpenFile(model.Parameters.MakeDir, outputPath, output.fileSuffix, modelName)
			if err != nil {
				println("unable to save "+name+": ", err)
			} else {
				rows := [][]string{output.columnNames}
				xColumnValue, yColumnValues, yLabels := output.values(model)
				rows = append(rows, append([]string{""}, yLabels...))
				for x := range xColumnValue {
					row := []string{strconv.FormatFloat(config.SIeV(xColumnValue[x], output.xUnit, model.Parameters.OutputUnits(), false), 'f', -1, 64)}
					for i := range yColumnValues[x] {
						row = append(row, strconv.FormatFloat(config.SIeV(yColumnValues[x][i], output.yUnit, model.Parameters.OutputUnits(), false), 'f', -1, 64))
					}
					rows = append(rows, row)
				}
				w := csv.NewWriter(file)
				w.WriteAll(rows)
				if model.Parameters.Verbose() {
					println(name + " saved")
				}
				if err := w.Error(); err != nil {
					log.Fatalln("error writing csv:", err)
				}
			}
		}
	}
}
