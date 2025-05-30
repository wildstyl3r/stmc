package utils

import (
	"encoding/csv"
	"os"
	"sort"

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
	scalarParams, err := OpenFile(true, path, subpath, clearName)
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
