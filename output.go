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
	saveFlag   *bool
	fileSuffix string
}

type SequentialDataItem struct {
	DataItem
	columnNames []string
	values      func(*DataExtractor) (args []float64, values [][]float64, labels []string)
	xUnit       []UnitElement
	yUnit       []UnitElement
}

type DataFlags struct {
	all         *bool
	sequentials map[string]SequentialDataItem
	outputPath  string
}

type DataExtractor struct {
	model       *Model
	cathodeFlux float64
	collisions  map[string][]float64
}

func newDataFlags() DataFlags {
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
					for x := range de.model.numCells {
						args = append(args, de.model.xStep*float64(x))
						values = append(values, []float64{de.model.VfromL(de.model.xStep * float64(x))})
					}
					return args, values, nil
				},
				xUnit: []UnitElement{{class: length, power: 1}},
				yUnit: []UnitElement{},
			},
			"Electric Field": {
				DataItem: DataItem{
					saveFlag:   flag.Bool("ef", false, "save Electric field"),
					fileSuffix: "Efield",
				},
				columnNames: []string{"x (cm)", "E (V/m)"},
				values: func(de *DataExtractor) (args []float64, values [][]float64, labels []string) {
					for x := range de.model.numCells {
						args = append(args, de.model.xStep*float64(x))
						values = append(values, []float64{de.model.EFieldFromL(de.model.xStep * float64(x))})
					}
					return args, values, nil
				},
				xUnit: []UnitElement{{class: length, power: 1}},
				yUnit: []UnitElement{{class: length, power: -1}},
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
				xUnit: []UnitElement{},
				yUnit: []UnitElement{{class: length, power: 1}},
			},
			"Collision counters": {
				DataItem: DataItem{
					saveFlag:   flag.Bool("cc", false, "save collision counters"),
					fileSuffix: "cc",
				},
				columnNames: []string{"x (cm)", "N_i(cm^{-1} Torr^{-1})"},
				values: func(de *DataExtractor) (args []float64, values [][]float64, labels []string) {
					for label := range de.collisions {
						labels = append(labels, label)
					}
					sort.Strings(labels)
					for x := range de.model.numCells {
						args = append(args, de.model.xStep*float64(x))
						var row []float64
						for _, label := range labels {
							row = append(row, de.collisions[label][x])
						}
						values = append(values, row)
					}
					return args, values, labels
				},
				xUnit: []UnitElement{{class: length, power: 1}},
				yUnit: []UnitElement{{class: length, power: -1}, {class: pressure, power: -1}},
			},
			"Normalized source term": {
				DataItem: DataItem{
					saveFlag:   flag.Bool("nst", true, "save normalized source term"),
					fileSuffix: "nst",
				},
				columnNames: []string{"x (cm)", "N_i(cm^{-1} Torr^{-1})"},
				values: func(de *DataExtractor) (args []float64, values [][]float64, labels []string) {
					labels = append(labels, string(lxgata.IONIZATION))
					sort.Strings(labels)
					for x := range de.model.numCells {
						args = append(args, de.model.xStep*float64(x))
						var row []float64
						row = append(row, de.collisions[string(lxgata.IONIZATION)][x])
						values = append(values, row)
					}
					return args, values, labels
				},
				xUnit: []UnitElement{{class: length, power: 1}},
				yUnit: []UnitElement{{class: length, power: -1}, {class: pressure, power: -1}},
			},
			"Energy loss due to ionizations": {
				DataItem: DataItem{
					saveFlag:   flag.Bool("li", false, "save ionization energy losses"),
					fileSuffix: "li",
				},
				columnNames: []string{"eV", "cm ^ -1"},
				values: func(de *DataExtractor) (args []float64, values [][]float64, labels []string) {
					for x := range de.model.numCells {
						args = append(args, de.model.xStep*float64(x))
						row := []float64{de.model.energyLossByProcess[string(lxgata.IONIZATION)][x]}
						values = append(values, row)
					}
					return args, values, []string{string(lxgata.IONIZATION)}
				},
				xUnit: []UnitElement{{class: length, power: 1}},
				yUnit: []UnitElement{{class: length, power: -1}},
			},
			"Out of energy for ionizations": {
				DataItem: DataItem{
					saveFlag:   flag.Bool("oe", false, "save ooe exit events count"),
					fileSuffix: "oe",
				},
				columnNames: []string{"x (cm)", "cm ^ -1"},
				values: func(de *DataExtractor) (args []float64, values [][]float64, labels []string) {
					sort.Strings(labels)
					for x := range de.model.numCells {
						args = append(args, de.model.xStep*float64(x))
						var row []float64
						row = append(row, float64((de.model.outOfEnergyAtCell[x]))/float64(de.model.parameters.NElectrons))
						values = append(values, row)
					}
					return args, values, nil
				},
				xUnit: []UnitElement{{class: length, power: 1}},
				yUnit: []UnitElement{{class: length, power: -1}},
			},
		},
	}
}

func newDataExtractor(model *Model) *DataExtractor {
	de := DataExtractor{
		model:       model,
		cathodeFlux: model.parameters.CathodeCurrentDensity / electronCharge, // [m^{-2} s^{-1}]
		collisions:  map[string][]float64{},
	}
	if de.model.parameters._dataFlags.outputPath != "" && de.model.parameters._dataFlags.outputPath[len(de.model.parameters._dataFlags.outputPath)-1] != '/' {
		de.model.parameters._dataFlags.outputPath += "/"
	}

	de.collisions = make(map[string][]float64)

	if model.parameters._verbose {
		collCounters := make(map[string]float64)

		for key, val := range model.collisionAtCell {
			if model.collisionAtCell[key] != nil {
				de.collisions[key] = make([]float64, de.model.numCells)
				for xIndex := range model.numCells {
					de.collisions[key][xIndex] = float64(val[xIndex]) / (float64(model.parameters.NElectrons) * model.xStep * model.parameters.Pressure)
					collCounters[key] += float64(model.collisionAtCell[key][xIndex])
				}
			}

		}
		for key := range collCounters {
			collCounters[key] /= float64(model.parameters.NElectrons)
		}
		fmt.Printf("Avg collisions per electron at distance %f: %v\n", model.parameters.CathodeFallLength, collCounters)
	}
	return &de
}

func openFile(makeDir bool, outputPath string, fileSuffix, modelName string) (*os.File, error) {
	if makeDir && fileSuffix != "" && fileSuffix != "." {
		os.MkdirAll(outputPath+fileSuffix, 0750)
		return os.Create(outputPath + fileSuffix + "/" + modelName + ".txt")
	} else {
		return os.Create(outputPath + modelName + "_" + fileSuffix + ".txt")
	}
}

func (de *DataExtractor) save(modelName string) {
	for name, output := range de.model.parameters._dataFlags.sequentials {
		if *output.saveFlag || *de.model.parameters._dataFlags.all {
			var file *os.File
			file, err := openFile(de.model.parameters.MakeDir, de.model.parameters._dataFlags.outputPath, output.fileSuffix, modelName)
			if err != nil {
				println("unable to save "+name+": ", err)
			} else {
				rows := [][]string{output.columnNames}
				xColumnValue, yColumnValues, yLabels := output.values(de)
				rows = append(rows, append([]string{""}, yLabels...))
				for x := range xColumnValue {
					row := []string{strconv.FormatFloat(SI(xColumnValue[x], output.xUnit, de.model.parameters._outputUnits, false), 'f', -1, 64)}
					for i := range yColumnValues[x] {
						row = append(row, strconv.FormatFloat(SI(yColumnValues[x][i], output.yUnit, de.model.parameters._outputUnits, false), 'f', -1, 64))
					}
					rows = append(rows, row)
				}
				w := csv.NewWriter(file)
				w.WriteAll(rows)
				if de.model.parameters._verbose {
					println(name + " saved")
				}
				if err := w.Error(); err != nil {
					log.Fatalln("error writing csv:", err)
				}
			}
		}
	}
}
