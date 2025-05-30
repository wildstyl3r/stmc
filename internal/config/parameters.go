package config

import (
	"fmt"
	"math"
	"os"
	"reflect"
	"slices"
	"strconv"
	"strings"

	"github.com/BurntSushi/toml"
	"github.com/wildstyl3r/lxgata"
	"github.com/wildstyl3r/stmc/internal/constants"
	"github.com/wildstyl3r/stmc/internal/utils"
)

type Config struct {
	OutputDir string
	Models    map[string]ModelParameters
	ModelParameters
	CVC                       string
	isDefinedMap              map[string]struct{}
	CurrentDividedByPressure2 bool
	CurrentDividedByArea      bool
	AddPotential              float64

	InputUnits  []string
	OutputUnits []string
}

func (c *Config) isDefined(path []string, meta *toml.MetaData) bool {
	if _, sureDefined := c.isDefinedMap[strings.Join(path, "#")]; sureDefined {
		return true
	} else {
		return meta.IsDefined(path...)
	}
}

func LoadConfig(configFileName string) (Config, toml.MetaData) {
	var config Config
	config.isDefinedMap = map[string]struct{}{}
	meta, err := toml.DecodeFile(configFileName+".toml", &config)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}

	var unitsConflict []string
	config.InputUnits, unitsConflict = checkUnits(config.InputUnits)
	if len(unitsConflict) > 0 {
		fmt.Printf("found input unit conflict: %v\n", unitsConflict)
		os.Exit(0)
	}
	if len(config.OutputUnits) == 0 {
		config.OutputUnits = config.InputUnits
	}
	config.OutputUnits, unitsConflict = checkUnits(config.OutputUnits)
	if len(unitsConflict) > 0 {
		fmt.Printf("found output unit conflict: %v\n", unitsConflict)
		os.Exit(0)
	}

	if len(config.CVC) > 0 {
		if len(config.Models) > 0 {
			fmt.Printf("Simultaneous CVC file listing and direct model specification not supported\n")
			os.Exit(0)
		}
		cvc, err := utils.ReadFloatPairs(config.CVC)
		if err != nil {
			fmt.Printf("CVC file reading error: %v\n", err)
			os.Exit(0)
		}
		filename := utils.GetFilename(config.CVC)
		config.Models = make(map[string]ModelParameters, len(cvc))
		for line := range cvc {
			modelName := filename + "_l" + strconv.Itoa(line+1)
			if config.CurrentDividedByPressure2 && !slices.Contains(valueUnits["CathodeCurrent"], UnitElement{Class: Pressure, Power: -2}) {
				valueUnits["CathodeCurrent"] = append(valueUnits["CathodeCurrent"], UnitElement{Class: Pressure, Power: -2})
			}
			if config.CurrentDividedByArea && !slices.Contains(valueUnits["CathodeCurrent"], UnitElement{Class: Length, Power: -2}) {
				valueUnits["CathodeCurrent"] = append(valueUnits["CathodeCurrent"], UnitElement{Class: Length, Power: -2})
			}
			config.Models[modelName] = ModelParameters{
				CathodeCurrent:       cvc[line][0],
				CathodeFallPotential: cvc[line][1],
			}
			config.isDefinedMap[strings.Join([]string{"Models", modelName, "CathodeCurrent"}, "#")] = struct{}{}
			config.isDefinedMap[strings.Join([]string{"Models", modelName, "CathodeFallPotential"}, "#")] = struct{}{}
		}
	} else {
		if len(config.Models) == 0 {
			fmt.Println("No models provided")
			os.Exit(0)
		}
	}

	return config, meta
}

type ModelParameters struct {
	CrossSections         string
	Species               string
	GapLength             float64 // [сm]
	PressureGapLength     float64 // [Pa*сm]
	CathodeFallLength     float64 // [сm]
	CathodeFallPotential  float64 // [V]
	CathodeCurrentDensity float64 // [mkA cm^{-2}] // ==>[A m^{-2}] == [C s^{-1} m^{-2}]
	CathodeCurrent        float64
	ConstEField           float64 // [V / m]
	Temperature           float64 // [K]
	Pressure              float64 // [Pa]
	// DAmbipolar              float64
	// AmbipolarCharacterScale float64
	CathodeRadius float64 // [cm]

	DeltaE                     float64 // [eV]
	NElectrons                 int
	MakeDir                    bool
	ParallelPlaneHollowCathode bool
	Volumetric                 bool
	CalculateStdError          bool

	CalculateCathodeFallLength bool
	CathodeFallLengthPrecision float64

	CountNulls bool

	GasDensity     float64
	_crossSections *lxgata.Collisions
	_outputUnits   []string
	_verbose       bool
	_threads       int
}

func (p *ModelParameters) CrossSectionsData() *lxgata.Collisions {
	return p._crossSections
}

func (p *ModelParameters) SetCrossSectionsData(cd *lxgata.Collisions) {
	p._crossSections = cd
}

func (p *ModelParameters) OutputUnits() []string {
	return p._outputUnits
}

func (p *ModelParameters) SetOutputUnits(u []string) {
	p._outputUnits = u
}

func (p *ModelParameters) Verbose() bool {
	return p._verbose
}

func (p *ModelParameters) SetVerbosity(verbose bool) {
	p._verbose = verbose
}

func (p *ModelParameters) Threads() int {
	return p._threads
}

func (p *ModelParameters) SetThreads(threads int) {
	p._threads = threads
}

var defaultValues = map[string]any{ // in SI
	"CathodeCurrentDensity":      0.0284,         //[A / m^2]
	"Pressure":                   101325. / 760., //[Pa]
	"ConstEField":                -100.,          //[V/m]
	"Temperature":                300.,           //[K]
	"CathodeFallLengthPrecision": 1e-6,           //[m]
	"DeltaE":                     0.1,
	"NElectrons":                 1000,
	"MakeDir":                    true,
	"ParallelPlaneHollowCathode": false,
	"CalculateCathodeFallLength": false,
	"Volumetric":                 false,
	"CountNulls":                 false,
}

var defaultUnits = []string{"mkA", "cm", "Torr"}

var fieldsXor = map[string][]string{
	"CalculateCathodeFallLength": {"CathodeFallLength"},
	"CathodeFallLength":          {"CalculateCathodeFallLength"},
	"PressureGapLength":          {"Pressure"},
	"Pressure":                   {"PressureGapLength"},
}

var fieldsAnd = map[string][]string{
	"Volumetric":                 {"CathodeRadius"},
	"CathodeCurrent":             {"CathodeRadius"},
	"CalculateCathodeFallLength": {"Species"},
	"PressureGapLength":          {"GapLength"},
}
var fieldsDerivable map[string][]string = map[string][]string{
	"CathodeCurrent":    {"CathodeCurrentDensity"},
	"PressureGapLength": {"Pressure"},
}

var valueUnits = map[string][]UnitElement{
	"GapLength": {
		{Class: Length, Power: 1},
	},
	"CathodeFallLength": {
		{Class: Length, Power: 1},
	},
	"CathodeCurrentDensity": {
		{Class: Current, Power: 1},
		{Class: Length, Power: -2},
	},
	"CathodeCurrent": {
		{Class: Current, Power: 1},
	},
	"ConstEField": {
		{Class: Length, Power: 1},
	},
	"PressureGapLength": {
		{Class: Pressure, Power: 1},
		{Class: Length, Power: 1},
	},
	"Pressure": {
		{Class: Pressure, Power: 1},
	},
	"CathodeRadius": {
		{Class: Length, Power: 1},
	},
}

var calculableFields = map[string]func(
	*ModelParameters,
	[]string,
) []string{
	"CathodeCurrent": func(mp *ModelParameters, definedFields []string) []string {
		if slices.Contains(definedFields, "CathodeRadius") || slices.Contains(valueUnits["CathodeCurrent"], UnitElement{Class: Length, Power: -2}) {
			temporary := mp.CathodeCurrent
			if slices.Contains(valueUnits["CathodeCurrent"], UnitElement{Class: Pressure, Power: -2}) {
				if slices.Contains(definedFields, "Pressure") {
					temporary *= mp.Pressure * mp.Pressure
				} else {
					return nil
				}
			}
			if !slices.Contains(valueUnits["CathodeCurrent"], UnitElement{Class: Length, Power: -2}) {
				area := math.Pi * mp.CathodeRadius * mp.CathodeRadius
				temporary /= area
			}
			mp.CathodeCurrentDensity = temporary
			return []string{"CathodeCurrentDensity"}
		}
		fmt.Printf("field 'CathodeRadius' not found: required by CathodeCurrentDensity calculation from CathodeCurrent\n")
		return nil
	},
	"PressureGapLength": func(mp *ModelParameters, definedFields []string) []string {
		if slices.Contains(definedFields, "GapLength") {
			mp.Pressure = mp.PressureGapLength / mp.GapLength
			return []string{"Pressure"}
		}
		fmt.Printf("field 'GapLength' not found: required by Pressure calculation from PressureGapLength\n")
		return nil
	},
}

func (modelConfig *ModelParameters) toSI(parameterNames, units []string) {
	modelConfigReflect := reflect.ValueOf(modelConfig).Elem()
	for name := range parameterNames {
		if modelConfigReflect.FieldByName(parameterNames[name]).CanFloat() {
			value := modelConfigReflect.FieldByName(parameterNames[name]).Float()
			value = SI(value, valueUnits[parameterNames[name]], units, true)
			modelConfigReflect.FieldByName(parameterNames[name]).SetFloat(value)
		}
	}
}

func (modelConfig *ModelParameters) checkFieldProblems(path []string, meta *toml.MetaData, globalConfig *Config) (ambiguities [][]string, missingDeps []string) {
	modelConfigReflect := reflect.ValueOf(modelConfig).Elem()
	for field := range fieldsXor {
		if globalConfig.isDefined(append(path, field), meta) {
			if modelConfigReflect.FieldByName(field).Kind() == reflect.Bool && !modelConfigReflect.FieldByName(field).Bool() {
				continue
			}
			var foundAlternatives []string
			for alternative := range fieldsXor[field] {
				if globalConfig.isDefined(append(path, fieldsXor[field][alternative]), meta) {
					foundAlternatives = append(foundAlternatives, fieldsXor[field][alternative])
				}
			}

			if len(foundAlternatives) > 0 {
				ambiguities = append(ambiguities, append([]string{field}, foundAlternatives...))
			}
		}
	}

	for field := range fieldsAnd {
		if globalConfig.isDefined(append(path, field), meta) {
			for requirement := range fieldsAnd[field] {
				if !globalConfig.isDefined(append(path, fieldsAnd[field][requirement]), meta) {
					missingDeps = append(missingDeps, fieldsAnd[field][requirement])
				}
			}
		}
	}
	return
}

/*
the algorithm:
0. preload into global and local
1. check problems in global
2. check problems in local
3. check combined
4. for local make list of exclusions from possible global & default
5. load missing from global
6. load missing from defaults
5. calculate calculables
6. check final missing
7. return success status

field value priority:
1. local
2. local-calculable
3. global
4. global-calculable
5. default
*/

func (modelConfig *ModelParameters) CheckAndUnify(modelName string, config *Config, meta *toml.MetaData) bool {
	globalAmbiguities, globalMissingDeps := config.checkFieldProblems([]string{}, meta, config)
	localAmbiguities, localMissingDeps := modelConfig.checkFieldProblems([]string{"Models", modelName}, meta, config)
	if len(globalAmbiguities) > 0 {
		fmt.Printf("unable to load config: found global ambiguities \n%v\n", globalAmbiguities)
		return false
	}
	if len(localAmbiguities) > 0 {
		fmt.Printf("unable to load config: found model ambiguities \n%v\n", localAmbiguities)
		return false
	}
	var missingIntersection []string
	for i := range globalMissingDeps {
		for j := range localMissingDeps {
			if globalMissingDeps[i] == localMissingDeps[j] {
				missingIntersection = append(missingIntersection, globalMissingDeps[i])
			}
		}
	}
	if len(missingIntersection) > 0 {
		fmt.Printf("unable to load config: required dependent fields not found \n%v\n", missingIntersection)
		return false
	}

	var discoveredParameters []string

	var excludeFromLoadingDefaultOrOuter map[string]struct{} = make(map[string]struct{})
	modelConfigReflect := reflect.ValueOf(&modelConfig).Elem()
	modelConfigType := modelConfigReflect.Elem().Type()
	for i := range modelConfigReflect.Elem().NumField() {
		fieldName := modelConfigType.Field(i).Name
		if config.isDefined([]string{"Models", modelName, fieldName}, meta) {
			discoveredParameters = append(discoveredParameters, fieldName)
			if xlist, some := fieldsXor[fieldName]; some {
				for x := range xlist {
					excludeFromLoadingDefaultOrOuter[xlist[x]] = struct{}{}
				}
			}
			if xlist, some := fieldsDerivable[fieldName]; some {
				for x := range xlist {
					excludeFromLoadingDefaultOrOuter[xlist[x]] = struct{}{}
				}
			}
		}
	}

	globalConfigReflect := reflect.ValueOf(&config).Elem()
	globalConfigType := globalConfigReflect.Elem().Type()
	for i := range globalConfigReflect.Elem().NumField() { // dive into embedded ModelParameters
		if globalConfigType.Field(i).Anonymous && globalConfigType.Field(i).Type.Kind() == reflect.Struct {
			globalConfigType = globalConfigReflect.Elem().Field(i).Type()
			globalConfigReflect = globalConfigReflect.Elem().Field(i)
			break
		}
	}

	for i := range globalConfigReflect.NumField() {
		fieldName := globalConfigType.Field(i).Name
		if _, some := excludeFromLoadingDefaultOrOuter[fieldName]; !some && !config.isDefined([]string{"Models", modelName, fieldName}, meta) && meta.IsDefined(fieldName) {
			modelConfigReflect.Elem().FieldByName(fieldName).Set(globalConfigReflect.Field(i))
			discoveredParameters = append(discoveredParameters, fieldName)
			excludeFromLoadingDefaultOrOuter[fieldName] = struct{}{}
			for xAlternative := range fieldsXor[fieldName] {
				excludeFromLoadingDefaultOrOuter[fieldsXor[fieldName][xAlternative]] = struct{}{}
			}
		}
	}

	modelConfig.toSI(discoveredParameters, config.InputUnits)

	for fieldName := range defaultValues {
		if _, x := excludeFromLoadingDefaultOrOuter[fieldName]; !x && !slices.Contains(discoveredParameters, fieldName) {
			modelConfigReflect.Elem().FieldByName(fieldName).Set(reflect.ValueOf(defaultValues[fieldName]))
			discoveredParameters = append(discoveredParameters, fieldName)
		}
	}

	var enabledParameters []string
	for fieldName := range discoveredParameters {
		field := modelConfigReflect.Elem().FieldByName(discoveredParameters[fieldName])
		if field.Kind() != reflect.Bool || field.Bool() {
			enabledParameters = append(enabledParameters, discoveredParameters[fieldName])
		}
	}

	calculatedAnything := true
	for calculatedAnything {
		calculatedAnything = false
		for initialFieldName := range calculableFields {
			if slices.Contains(enabledParameters, initialFieldName) {
				calculated := calculableFields[initialFieldName](modelConfig, enabledParameters)
				if len(calculated) != 0 {
					calculatedAnything = true
					enabledParameters = append(enabledParameters, calculated...)
					enabledParameters = slices.DeleteFunc(enabledParameters, func(elem string) bool {
						return elem == initialFieldName
					})
				}
			}
		}
	}
	for initialFieldName := range calculableFields {
		if slices.Contains(enabledParameters, initialFieldName) {
			return false
		}
	}

	allGood := true

	for i := range enabledParameters {
		for requirement := range fieldsAnd[enabledParameters[i]] {
			if !slices.Contains(enabledParameters, fieldsAnd[enabledParameters[i]][requirement]) {
				fmt.Printf("for parameter %s requirement %s not found\n", enabledParameters[i], fieldsAnd[enabledParameters[i]][requirement])
				allGood = false
			}
		}
		for conflict := range fieldsXor[enabledParameters[i]] {
			if slices.Contains(enabledParameters, fieldsXor[enabledParameters[i]][conflict]) {
				fmt.Printf("for parameter %s found conflicting parameter: %s\n", enabledParameters[i], fieldsXor[enabledParameters[i]][conflict])
				allGood = false
			}
		}
	}

	modelConfig.CathodeFallPotential += config.AddPotential
	modelConfig.GasDensity = modelConfig.Pressure / (constants.KBolzmann * modelConfig.Temperature)
	var conflict []string
	units, conflict := checkUnits(config.OutputUnits)
	if len(conflict) > 0 {
		fmt.Printf("found output unit conflict: %v\n Data will be saved in input units", conflict)
		modelConfig._outputUnits = config.InputUnits
	} else {
		modelConfig._outputUnits = units
	}

	return allGood
}
