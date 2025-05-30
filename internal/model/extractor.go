package model

import (
	"encoding/csv"
	"fmt"
	"log"
	"math"
	"os"
	"strconv"

	"github.com/wildstyl3r/lxgata"
	"github.com/wildstyl3r/stmc/internal/config"
	"github.com/wildstyl3r/stmc/internal/constants"
	"github.com/wildstyl3r/stmc/internal/utils"
)

type DataExtractor struct {
	model           *Model
	cathodeFlux     float64
	collisions      map[lxgata.CollisionType][]float64
	collisionsError map[lxgata.CollisionType][]float64
}

func NewDataExtractor(model *Model) *DataExtractor {
	de := DataExtractor{
		model:           model,
		cathodeFlux:     model.Parameters.CathodeCurrentDensity / constants.ElectronCharge, // [m^{-2} s^{-1}]
		collisions:      map[lxgata.CollisionType][]float64{},
		collisionsError: map[lxgata.CollisionType][]float64{},
	}

	// de.collisions = make(map[lxgata.CollisionType][]float64)

	if model.Parameters.Verbose() {
		collCounters := make(map[lxgata.CollisionType]float64)

		for key, val := range model.CollisionAtCell {
			if model.CollisionAtCell[key] != nil {
				de.collisions[key] = make([]float64, de.model.NumCells)
				de.collisionsError[key] = make([]float64, de.model.NumCells)
				for xIndex := range model.NumCells {
					de.collisions[key][xIndex] = float64(utils.SumSlice(val[xIndex])) / (float64(model.Parameters.NElectrons) * model.XStep * model.Parameters.Pressure)
					de.collisionsError[key][xIndex] = constants.Quantile95 * math.Sqrt(float64(utils.Variance(val[xIndex], true))) / (math.Sqrt(float64(model.Parameters.NElectrons)) * model.XStep * model.Parameters.Pressure)
					// 1.96 is double-sided quantile for 95% confidence
					collCounters[key] += float64(utils.SumSlice(model.CollisionAtCell[key][xIndex]))
				}
			}

		}
		for key := range collCounters {
			collCounters[key] /= float64(model.Parameters.NElectrons)
		}
		fmt.Printf("Avg collisions per electron at distance %f: %v\n", model.Parameters.CathodeFallLength, collCounters)
	}
	return &de
}

func (de *DataExtractor) Save(modelName string, df DataFlags) {
	for name, output := range df.sequentials {
		if *output.saveFlag || *df.all {
			var file *os.File
			file, err := utils.OpenFile(de.model.Parameters.MakeDir, df.outputPath, output.fileSuffix, modelName)
			if err != nil {
				println("unable to save "+name+": ", err)
			} else {
				rows := [][]string{output.columnNames}
				xColumnValue, yColumnValues, yLabels := output.values(de)
				rows = append(rows, append([]string{""}, yLabels...))
				for x := range xColumnValue {
					row := []string{strconv.FormatFloat(config.SI(xColumnValue[x], output.xUnit, de.model.Parameters.OutputUnits(), false), 'f', -1, 64)}
					for i := range yColumnValues[x] {
						row = append(row, strconv.FormatFloat(config.SI(yColumnValues[x][i], output.yUnit, de.model.Parameters.OutputUnits(), false), 'f', -1, 64))
					}
					rows = append(rows, row)
				}
				w := csv.NewWriter(file)
				w.WriteAll(rows)
				if de.model.Parameters.Verbose() {
					println(name + " saved")
				}
				if err := w.Error(); err != nil {
					log.Fatalln("error writing csv:", err)
				}
			}
		}
	}
}
