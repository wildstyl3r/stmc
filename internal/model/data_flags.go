package model

import (
	"flag"
	"sort"

	"github.com/wildstyl3r/lxgata"
	"github.com/wildstyl3r/stmc/internal/config"
)

type DataItem struct {
	saveFlag   *bool
	fileSuffix string
}

type SequentialDataItem struct {
	DataItem
	columnNames []string
	values      func(*DataExtractor) (args []float64, values [][]float64, labels []string)
	xUnit       []config.UnitElement
	yUnit       []config.UnitElement
}

type DataFlags struct {
	all         *bool
	sequentials map[string]SequentialDataItem
	outputPath  string
}

func NewDataFlags() DataFlags {
	return DataFlags{
		all: flag.Bool("all", false, "save every available metric"),
		sequentials: map[string]SequentialDataItem{
			"Potential": {
				DataItem: DataItem{
					saveFlag:   flag.Bool("p", false, "save potential"),
					fileSuffix: "V",
				},
				columnNames: []string{"x (cm)", "g (V)"},
				values: func(de *DataExtractor) (args []float64, values [][]float64, labels []string) {
					for x := range de.model.NumCells {
						args = append(args, de.model.XStep*float64(x))
						values = append(values, []float64{de.model.VfromL(de.model.XStep * float64(x))})
					}
					return args, values, nil
				},
				xUnit: []config.UnitElement{{Class: config.Length, Power: 1}},
				yUnit: []config.UnitElement{},
			},
			"Electric Field": {
				DataItem: DataItem{
					saveFlag:   flag.Bool("ef", false, "save Electric field"),
					fileSuffix: "Efield",
				},
				columnNames: []string{"x (cm)", "E (V/m)"},
				values: func(de *DataExtractor) (args []float64, values [][]float64, labels []string) {
					for x := range de.model.NumCells {
						args = append(args, de.model.XStep*float64(x))
						values = append(values, []float64{de.model.EFieldFromL(de.model.XStep * float64(x))})
					}
					return args, values, nil
				},
				xUnit: []config.UnitElement{{Class: config.Length, Power: 1}},
				yUnit: []config.UnitElement{{Class: config.Length, Power: -1}},
			},
			"LfromV": {
				DataItem: DataItem{
					saveFlag:   flag.Bool("lv", false, "save x from v"),
					fileSuffix: "lv",
				},
				columnNames: []string{"g (V)", "x (cm)"},
				values: func(de *DataExtractor) (args []float64, values [][]float64, labels []string) {
					for i := range 100000 {
						v := -(de.model.Va + de.model.Vc) * float64(100000-i) / 100000.
						args = append(args, v)
						values = append(values, []float64{de.model.LfromV(v)})
					}
					return args, values, nil
				},
				xUnit: []config.UnitElement{},
				yUnit: []config.UnitElement{{Class: config.Length, Power: 1}},
			},
			"Collision counters": {
				DataItem: DataItem{
					saveFlag:   flag.Bool("cc", false, "save collision counters"),
					fileSuffix: "cc",
				},
				columnNames: []string{"x (cm)", "N_i(cm^{-1} Torr^{-1})", "Standard error"},
				values: func(de *DataExtractor) (args []float64, values [][]float64, labels []string) {
					for label := range de.collisions {
						labels = append(labels, string(label))
						if de.model.Parameters.CalculateStdError {
							labels = append(labels, string(label)+"_conf_interval")
						}
					}
					sort.Strings(labels)
					for x := range de.model.NumCells {
						args = append(args, de.model.XStep*float64(x))
						var row []float64
						for _, label := range labels {
							row = append(row, de.collisions[lxgata.CollisionType(label)][x])
							if de.model.Parameters.CalculateStdError {
								row = append(row, de.collisionsError[lxgata.CollisionType(label)][x])
							}
						}
						values = append(values, row)
					}
					return args, values, labels
				},
				xUnit: []config.UnitElement{{Class: config.Length, Power: 1}},
				yUnit: []config.UnitElement{{Class: config.Length, Power: -1}, {Class: config.Pressure, Power: -1}},
			},
			"Normalized source term": {
				DataItem: DataItem{
					saveFlag:   flag.Bool("nst", true, "save normalized source term"),
					fileSuffix: "nst",
				},
				columnNames: []string{"x (cm)", "N_i(cm^{-1} Torr^{-1})", "Standard error"},
				values: func(de *DataExtractor) (args []float64, values [][]float64, labels []string) {
					labels = append(labels, string(lxgata.IONIZATION))
					if de.model.Parameters.CalculateStdError {
						labels = append(labels, string(lxgata.IONIZATION)+"_conf_interval")
					}
					sort.Strings(labels)
					for x := range de.model.NumCells {
						args = append(args, de.model.XStep*float64(x))
						var row []float64
						row = append(row, de.collisions[lxgata.IONIZATION][x])
						if de.model.Parameters.CalculateStdError {
							row = append(row, de.collisionsError[lxgata.IONIZATION][x])
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
					saveFlag:   flag.Bool("li", false, "save ionization energy losses"),
					fileSuffix: "li",
				},
				columnNames: []string{"eV", "cm ^ -1"},
				values: func(de *DataExtractor) (args []float64, values [][]float64, labels []string) {
					for x := range de.model.NumCells {
						args = append(args, de.model.XStep*float64(x))
						row := []float64{de.model.EnergyLossByProcess[lxgata.IONIZATION][x]}
						values = append(values, row)
					}
					return args, values, []string{string(lxgata.IONIZATION)}
				},
				xUnit: []config.UnitElement{{Class: config.Length, Power: 1}},
				yUnit: []config.UnitElement{{Class: config.Length, Power: -1}},
			},
			"Out of energy for ionizations": {
				DataItem: DataItem{
					saveFlag:   flag.Bool("oe", false, "save ooe exit events count"),
					fileSuffix: "oe",
				},
				columnNames: []string{"x (cm)", "cm ^ -1"},
				values: func(de *DataExtractor) (args []float64, values [][]float64, labels []string) {
					for label := range de.collisions {
						labels = append(labels, string(label))
						if de.model.Parameters.CalculateStdError {
							labels = append(labels, string(label)+"_conf_interval")
						}
					}
					sort.Strings(labels)
					for x := range de.model.NumCells {
						args = append(args, de.model.XStep*float64(x))
						var row []float64
						row = append(row, float64((de.model.OutOfEnergyAtCell[x]))/float64(de.model.Parameters.NElectrons))
						values = append(values, row)
					}
					return args, values, nil
				},
				xUnit: []config.UnitElement{{Class: config.Length, Power: 1}},
				yUnit: []config.UnitElement{{Class: config.Length, Power: -1}},
			},
			"Plasma density": {
				DataItem: DataItem{
					saveFlag:   flag.Bool("n", true, "save plasma density"),
					fileSuffix: "n",
				},
				columnNames: []string{"x (cm)", "cm ^ -3"},
				values: func(de *DataExtractor) (args []float64, values [][]float64, labels []string) {
					density := GlowDischargeDensity(de)
					for x := range de.model.NumCells {
						args = append(args, de.model.XStep*float64(x))
						values = append(values, []float64{density[x]})
					}
					return args, values, []string{"Plasma density n(x)"}
				},
				xUnit: []config.UnitElement{{Class: config.Length, Power: 1}},
				yUnit: []config.UnitElement{{Class: config.Length, Power: -3}},
			},
		},
	}
}

func (df *DataFlags) SetOutputPath(path string) {

	if path != "" && path[len(path)-1] != '/' {
		df.outputPath = path + "/"
	} else {
		df.outputPath = path
	}
}

func (df *DataFlags) GetOutputPath() string {
	return df.outputPath
}
