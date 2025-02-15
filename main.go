package main

import (
	"encoding/csv"
	"flag"
	"fmt"
	"log"
	"math"
	"os"
	"runtime"
	"runtime/pprof"
	"strconv"
	"strings"
	"time"

	"github.com/wildstyl3r/lxgata"
)

func main() {
	cpuprofile := flag.String("cpuprofile", "", "write cpu profile to file")
	threads := flag.Int("j", runtime.NumCPU(), "threads to run")
	dataExtractorFlags := newDataFlags()
	configFileNamePointer := flag.String("input", "inputs/mini_st/miniHe", "model configuration in toml format") //"inputs/val/BM_He", "model configuration in toml format")
	flag.Parse()

	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	var config, meta = loadConfig(strings.TrimSuffix(*configFileNamePointer, ".toml"))

	if config.OutputDir != "" && config.OutputDir != "." {
		os.MkdirAll(config.OutputDir, 0750)
		dataExtractorFlags.outputPath += config.OutputDir + "/"
	}

	var scalarW *csv.Writer
	println(*configFileNamePointer, dataExtractorFlags.outputPath)
	nameSegments := strings.Split(strings.TrimSuffix(*configFileNamePointer, ".toml"), "/")
	clearName := nameSegments[len(nameSegments)-1]
	scalarParams, err := openFile(true, dataExtractorFlags.outputPath, "calculated_params", clearName)
	if err != nil {
		println("unable to save dc and secondary emission coefficient: ", err.Error())
		os.Exit(1)
	} else {
		scalarW = csv.NewWriter(scalarParams)
		scalarW.WriteAll([][]string{
			{"E/N", "integrated secondary emission coefficient", "analytic secondary emission coefficient", "final gamma loss", "sheath length"},
		})
		scalarW.Flush()
	}

	// var lossW *csv.Writer
	// lossFile, err := openFile(true, dataExtractorFlags.outputPath, "loss_debug", clearName)
	// if err != nil {
	// 	println("unable to save loss debug data: ", err.Error())
	// 	os.Exit(1)
	// } else {
	// 	lossW = csv.NewWriter(lossFile)
	// 	lossW.WriteAll([][]string{
	// 		{"x, cm", "loss"},
	// 	})
	// 	lossW.Flush()
	// }

	crossSections, err := lxgata.LoadCrossSections(config.CrossSections)
	if err != nil {
		panic(fmt.Errorf("invalid cross section file: %w", err))
	}
	for i := range crossSections {
		crossSections[i].Expand(1e-4)
	}

	for modelName, parameters := range config.Models {
		startTime := time.Now()
		fmt.Printf("Current time: %s\n", startTime.UTC().Format(time.UnixDate))

		runtime.GC()
		fmt.Println("\n" + modelName)
		parameters.checkMetrizeSetDefaults(modelName, &config, &meta)

		if parameters.CalculateCathodeFallLength {
			var gamma_loss, g_integral, g_analytic float64
			fmt.Print("calculating dc iteratively:\n")
			iteration := 1
			itp := &iteration
			density := parameters.Pressure / (k * parameters.Temperature)
			// for i := 1; i < int(config.GapLength/0.0005); i++ {
			// 	dc := float64(i) * 0.0005
			// 	fmt.Printf("step %d\n", *itp)
			// 	*itp += 1
			// 	model := Model{
			// 		crossSections: crossSections,
			// 		dataFlags:     dataExtractorFlags,
			// 		Vc:            -math.Abs(parameters.CathodeFallPotential),
			// 		parameters:    parameters,
			// 		gasDensity:    density,
			// 		threads:       *threads,
			// 	}
			// 	model.parameters.CathodeFallLength = dc
			// 	model.init()
			// 	model.run()
			// 	dataExtractor := DataExtractor{
			// 		flags: dataExtractorFlags,
			// 	}
			// 	dataExtractor.extract(modelName, &model, &parameters)
			// 	dataExtractor.save(modelName, &parameters)

			// 	var gamma_loss float64
			// 	gamma_loss, g_integral, g_analytic = dataExtractor.gamma_loss()
			// 	fmt.Printf("d_c: %f.6\nsecondary emission coefficient\n\t integral: %f.6\n\t analytic:%f.6\n", dc, g_integral, g_analytic)
			// 	lossW.WriteAll([][]string{
			// 		{strconv.FormatFloat(dc, 'f', 10, 64), strconv.FormatFloat(gamma_loss, 'f', 10, 64)},
			// 	})
			// }
			dc := ternarySearchMax(func(dc float64) float64 {
				fmt.Printf("step %d\n", *itp)
				*itp += 1
				model := Model{
					crossSections: &crossSections,
					dataFlags:     dataExtractorFlags,
					Vc:            math.Abs(parameters.CathodeFallPotential),
					parameters:    parameters,
					gasDensity:    density,
					threads:       *threads,
				}
				model.parameters.CathodeFallLength = dc
				model.init()
				model.run()
				dataExtractor := DataExtractor{
					flags: dataExtractorFlags,
				}
				dataExtractor.extract(modelName, &model, &parameters)
				dataExtractor.save(modelName, &parameters)
				gamma_loss, g_integral, g_analytic = dataExtractor.gamma_loss()
				fmt.Printf("d_c: %6f\nsecondary emission coefficient\n\t integral: %6f\n\t analytic:%6f\n", dc, g_integral, g_analytic)
				return -gamma_loss
			}, 0, parameters.GapLength, parameters.CathodeFallLengthPrecision)
			scalarW.WriteAll([][]string{
				{strconv.FormatFloat(parameters.CathodeFallPotential/(dc*density*Townsend), 'f', 10, 64), strconv.FormatFloat(g_integral, 'f', 10, 64), strconv.FormatFloat(g_analytic, 'f', 10, 64), strconv.FormatFloat(gamma_loss, 'f', 10, 64), strconv.FormatFloat(dc, 'f', 10, 64)},
			})
			println("saved d_c and γ")
			scalarW.Flush()

		} else {
			model := Model{
				crossSections: &crossSections,
				dataFlags:     dataExtractorFlags,
				Vc:            math.Abs(parameters.CathodeFallPotential),
				parameters:    parameters,
				gasDensity:    parameters.Pressure / (k * parameters.Temperature),
				threads:       *threads,
			}
			model.init()
			model.run()
			dataExtractor := DataExtractor{
				flags: dataExtractorFlags,
			}
			dataExtractor.extract(modelName, &model, &parameters)

			dataExtractor.save(modelName, &parameters)
		}
		fmt.Printf("Elapsed time: %v\n", time.Since(startTime))
	}
}
