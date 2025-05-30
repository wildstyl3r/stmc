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
	"github.com/wildstyl3r/stmc/internal/config"
	"github.com/wildstyl3r/stmc/internal/constants"
	"github.com/wildstyl3r/stmc/internal/model"
	"github.com/wildstyl3r/stmc/internal/utils"
)

func main() {
	startTime := time.Now()
	cpuprofile := flag.String("cpuprofile", "", "write cpu profile to file")
	memprofile := flag.String("memprofile", "", "write memory profile to `file`")
	threads := flag.Int("j", runtime.NumCPU(), "threads to run")
	verbose := flag.Bool("v", true, "verbose")
	debug := flag.Bool("d", false, "debug")
	dataExtractorFlags := model.NewDataFlags()                                                                           //Donko/cvc_45mbarcm3Dup
	configFileNamePointer := flag.String("i", "inputs/donko2009/240Pa_cm_3D.toml", "model configuration in toml format") //"inputs/val/BM_He", "model configuration in toml format")
	rootFindingAlgorithm := flag.String("alg", "s", "root finder algorithm for gamma calculation ([s]tochastic approximation, [b]inary search, [t]ernary search)")
	flag.Parse()

	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}
	fmt.Printf("CONFIG: %s\n", *configFileNamePointer)
	var config, meta = config.LoadConfig(strings.TrimSuffix(*configFileNamePointer, ".toml"))

	if config.OutputDir != "" && config.OutputDir != "." {
		outputDir := strings.TrimSuffix(config.OutputDir, "/") + "/" + utils.GetFilename(*configFileNamePointer) + "/"
		fmt.Printf("OUTPUT DIR: %s\n", outputDir)
		os.MkdirAll(config.OutputDir, 0750)
		dataExtractorFlags.SetOutputPath(outputDir)
	} else {
		panic(fmt.Errorf("output folder not specified"))
	}

	crossSections, err := lxgata.LoadCrossSections(config.CrossSections)
	if err != nil {
		panic(fmt.Errorf("invalid cross section file: %w", err))
	}
	var gammaData utils.CSV

	for modelName, parameters := range config.Models {
		parameters.SetCrossSectionsData(&crossSections)
		parameters.SetVerbosity(*verbose)
		parameters.SetThreads(*threads)

		modelStartTime := time.Now()
		fmt.Printf("Current time: %s\n", modelStartTime.UTC().Format(time.UnixDate))

		runtime.GC()
		fmt.Println("\n" + modelName)
		if !parameters.CheckAndUnify(modelName, &config, &meta) {
			fmt.Printf("found a problem in the config for [%v], skipping\n", modelName)
			continue
		}

		if parameters.CalculateCathodeFallLength {
			var gammaLoss, gammaIntegral, gammaVariance, gammaCI, gammaAnalytic float64
			if *verbose {
				fmt.Print("calculating dc iteratively:\n")
			}
			iteration := 1
			itp := &iteration
			minDc, maxDc := model.EstimateCathodeFallLengthLimits(&parameters)
			if *debug {
				dcStep := 0.00001
				nSteps := int((maxDc - minDc) / dcStep)
				fmt.Printf("GapLen: %f, nSteps: %d\n", config.GapLength, nSteps)
				gamma := make([]float64, nSteps)
				intS := make([]float64, nSteps)
				debugData := [][]string{{"dc", "E/N", "integrated secondary emission coefficient", "analytic secondary emission coefficient", "gamma difference"}}
				for i := range gamma {
					dc := float64(i)*dcStep + minDc
					var lossType utils.LossType
					switch *rootFindingAlgorithm {
					case "t":
						lossType = utils.MSE
					case "b", "s":
						lossType = utils.Difference
					}
					gammaLoss, gamma[i], gammaAnalytic, _ = model.GammaCalculationStep(itp, dc, parameters, lossType)
					intS[i] = 1. / gamma[i]
					debugData = append(debugData, []string{
						strconv.FormatFloat(dc, 'f', 10, 64),
						strconv.FormatFloat(parameters.CathodeFallPotential/(dc*parameters.GasDensity*constants.Townsend), 'f', 10, 64),
						strconv.FormatFloat(gamma[i], 'f', 10, 64),
						strconv.FormatFloat(gammaAnalytic, 'f', 10, 64),
						strconv.FormatFloat(gammaLoss, 'f', 10, 64), //strconv.FormatFloat(gammaVariance, 'f', 10, 64),
					})
				}
				gammaMean, gammaVariance := utils.MeanAndVariance(gamma, true)
				intSMean, intSVariance := utils.MeanAndVariance(intS, true)
				fmt.Printf("[%s] gamma mean: %.9f, gamma variance: %.9f, integral S mean: %.9f, integral S variance: %.9f\n", modelName, gammaMean, gammaVariance, intSMean, intSVariance)
				debugFile, err := utils.OpenFile(true, dataExtractorFlags.GetOutputPath(), "debug", utils.GetFilename(*configFileNamePointer)+modelName)

				if err != nil {
					println("unable to save dc and secondary emission coefficient: ", err.Error())
					os.Exit(1)
				} else {
					debugW := csv.NewWriter(debugFile)
					debugW.WriteAll(debugData)
					debugW.Flush()
				}
			} else {
				var dc float64
				var dataExtractor *model.DataExtractor
				switch *rootFindingAlgorithm {
				case "s":
					fmt.Println("calculating gamma with stochastic approximation")
					var fLeft, fRight, initialDc float64 // := f(left+0.5*thetaPrecision), f(right-0.5*thetaPrecision)
					{
						var gILeft, gIRight float64
						fLeft, gILeft, _, _ = model.GammaCalculationStep(nil, minDc+parameters.CathodeFallLengthPrecision, parameters, utils.Difference)
						fRight, gIRight, _, _ = model.GammaCalculationStep(nil, maxDc-parameters.CathodeFallLengthPrecision, parameters, utils.Difference)

						meanGI := 0.5 * (gILeft + gIRight)
						initialDc = model.GetApproximateDcForGamma(meanGI, minDc, maxDc, parameters)
					}
					approxLossDerivative := (fRight - fLeft) / (maxDc - minDc - 2*parameters.CathodeFallLengthPrecision)

					gammaI := []float64{}
					dc = utils.StochasticApproximation(minDc, maxDc, initialDc, approxLossDerivative, parameters.CathodeFallLengthPrecision, constants.Quantile95, 10, func(dc float64) float64 {
						gammaLoss, gammaIntegral, gammaAnalytic, dataExtractor = model.GammaCalculationStep(itp, dc, parameters, utils.Difference)
						gammaI = append(gammaI, gammaIntegral)
						return gammaLoss
					})
					gammaIntegral, gammaVariance = utils.MeanAndVariance(gammaI[len(gammaI)-10:], true)
					gammaCI = 2 * math.Sqrt(gammaVariance) * constants.Quantile95 / math.Sqrt(float64(10))
					gammaLoss = gammaAnalytic - gammaIntegral

				case "b":
					fmt.Println("calculating gamma with naive bisection")
					dcLeft, dcRight := utils.BinarySearch(func(dc float64) bool {
						gammaLoss, gammaIntegral, gammaAnalytic, dataExtractor = model.GammaCalculationStep(itp, dc, parameters, utils.Difference)
						return gammaLoss > 0
					}, minDc, min(maxDc, parameters.GapLength), parameters.CathodeFallLengthPrecision) //0, parameters.GapLength, parameters.CathodeFallLengthPrecision)
					dc = 0.5 * (dcLeft + dcRight)

				case "t":
					fmt.Println("calculating gamma with naive ternary search")
					fmt.Printf("min dc: %f, max dc: %f, gap len: %f\n", minDc, maxDc, parameters.GapLength)
					dc = utils.TernarySearchMax(func(dc float64) float64 {
						gammaLoss, gammaIntegral, gammaAnalytic, dataExtractor = model.GammaCalculationStep(itp, dc, parameters, utils.MSE)
						return -gammaLoss
					}, minDc, min(maxDc, parameters.GapLength), parameters.CathodeFallLengthPrecision) //0, parameters.GapLength, parameters.CathodeFallLengthPrecision)
				}
				dataExtractor.Save(modelName, dataExtractorFlags)
				gammaData = append(gammaData,
					[]string{
						strconv.FormatFloat(parameters.CathodeFallPotential/(dc*parameters.GasDensity*constants.Townsend), 'f', 10, 64),
						strconv.FormatFloat(gammaIntegral, 'f', 10, 64),
						strconv.FormatFloat(gammaAnalytic, 'f', 10, 64),
						strconv.FormatFloat(gammaLoss, 'f', 10, 64),
						strconv.FormatFloat(dc, 'f', 10, 64),
						strconv.FormatFloat(gammaCI, 'f', 10, 64),
					})
				println("saved d_c and γ")
			}

		} else {
			m := model.NewModel(parameters.CathodeFallLength, parameters)
			m.Run()
			dataExtractor := model.NewDataExtractor(&m)

			dataExtractor.Save(modelName, dataExtractorFlags)
		}
		fmt.Printf("Elapsed time: %v\n", time.Since(modelStartTime))
	}
	if len(gammaData) > 0 {
		utils.WriteAsCSV(
			gammaData,
			dataExtractorFlags.GetOutputPath(), "gamma", *configFileNamePointer,
			[]string{"E/N", "integrated secondary emission coefficient", "analytic secondary emission coefficient", "final gamma loss", "sheath length", "integrated secondary emission coefficient_conf_interval"},
		)
	}
	fmt.Printf("Total elapsed time: %v\n\n", time.Since(startTime))
	if *memprofile != "" {
		f, err := os.Create(*memprofile)
		if err != nil {
			log.Fatal("could not create memory profile: ", err)
		}
		defer f.Close() // error handling omitted for example
		runtime.GC()    // get up-to-date statistics
		// Lookup("allocs") creates a profile similar to go test -memprofile.
		// Alternatively, use Lookup("heap") for a profile
		// that has inuse_space as the default index.
		if err := pprof.Lookup("allocs").WriteTo(f, 0); err != nil {
			log.Fatal("could not write memory profile: ", err)
		}
	}
}
