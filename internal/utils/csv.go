package utils

import (
	"encoding/csv"
	"fmt"
	"os"
	"sort"
	"strconv"

	"github.com/facette/natsort"
)

type CSV [][]string

func (data CSV) Less(i, j int) bool {
	return natsort.Compare(data[i][0], data[j][0])
}

func (data CSV) Len() int {
	return len(data)
}
func (data CSV) Swap(i, j int) {
	data[i], data[j] = data[j], data[i]
}

func WriteAsCSV(data CSV, path, subpath, filename string, columns []string) {
	var scalarW *csv.Writer
	println(filename, path)
	clearName := GetFilename(filename)
	scalarParams, err := OpenFile(true, path, subpath, clearName, TypeTXT)
	if err != nil {
		println("unable to save dc and secondary emission coefficient: ", err.Error())
		os.Exit(1)
	} else {
		scalarW = csv.NewWriter(scalarParams)
		scalarW.WriteAll([][]string{
			columns,
		})
		scalarW.Flush()
	}
	sort.Sort(data)
	scalarW.WriteAll(data)
	scalarW.Flush()
}

func ReadDataCSV(filename string) (metadata [][]string, columnNames []string, values [][]float64, err error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, nil, nil, fmt.Errorf("error opening file: %w", err)
	}
	defer file.Close()

	reader := csv.NewReader(file)

	rows, err := reader.ReadAll()
	if err != nil {
		return nil, nil, nil, fmt.Errorf("error reading file: %w", err)
	}

	getData := func(row []string) (dataRow []float64) {
		dataRow = make([]float64, len(row))
		for i := range row {
			value, err := strconv.ParseFloat(row[i], 64)
			if err != nil {
				return nil
			}
			dataRow[i] = value
		}
		return
	}

	dataWidth := 0
	for i := range rows {
		data := getData(rows[i])
		if data != nil {
			if i == 0 {
				return nil, nil, nil, fmt.Errorf("no column names provided, unable to inference column meaning")
			}
			if dataWidth == 0 {
				dataWidth = len(data)
				if i > 1 {
					metadata = rows[:i-2]
				}
				columnNames = rows[i-1]
			} else if dataWidth != len(data) {
				return nil, nil, nil, fmt.Errorf("wrong number of columns: %d vs previous %d at line %d", len(data), dataWidth, i)
			}
			values = append(values, data)
		} else if dataWidth != 0 {
			return nil, nil, nil, fmt.Errorf("unexpected non-numerics beyond data section at line %d", i)
		}
	}
	return
}
