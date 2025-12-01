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
	"github.com/gocarina/gocsv"
	"github.com/wildstyl3r/lxgata"
	"github.com/wildstyl3r/stmc/internal/constants"
	"github.com/wildstyl3r/stmc/internal/utils"
)

type IonizationEnergySharing int

const (
	Zero IonizationEnergySharing = iota
	UniformRandom
	Equal
	Opal
)

type IonizationScattering int

const (
	Boeuf IonizationScattering = iota
	Isotropic
	Coulomb
	CoulombM
	Born
	BornFullLoss
	Mixed
	MixedFullLoss
)

type UnitConfig struct {
	L string
	T string
	P string
	I string
	E string
}

func (uc *UnitConfig) GetForClass(class UnitClass) string {
	switch class {
	case Current:
		return uc.I
	case Energy:
		return uc.E
	case Length:
		return uc.L
	case Pressure:
		return uc.P
	case Time:
		return uc.T
	case Voltage:
		return "V"
	default:
		return ""
	}
}

type Config struct {
	OutputDir string
	Models    map[string]ModelParameters
	Options   map[string]map[string][]string
	ModelParameters
	ModelTable   string
	isDefinedMap map[string]struct{}
	AddPotential float64

	InputUnits  UnitConfig
	OutputUnits UnitConfig
}

func (c *Config) isDefined(path []string, meta *toml.MetaData) bool {
	if len(path) == 3 { // the only possible values are 1 for global and 3 for per-model definition
		// model might be renamed by Options
		prototypePath := []string{path[0], c.Models[path[1]]._prototypeName, path[2]}
		if _, sureDefined := c.isDefinedMap[strings.Join(prototypePath, "#")]; sureDefined {
			return true
		}
		if meta.IsDefined(prototypePath...) {
			return true
		}
	}
	if _, sureDefined := c.isDefinedMap[strings.Join(path, "#")]; sureDefined {
		return true
	} else {
		return meta.IsDefined(path...)
	}
}

func LoadConfig(flags Flags) (Config, toml.MetaData) {
	var config Config
	config.isDefinedMap = map[string]struct{}{}
	meta, err := toml.DecodeFile(strings.TrimSuffix(*flags.ConfigFileNamePointer, ".toml")+".toml", &config)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}

	var unknownUnits []string
	config.InputUnits, unknownUnits = checkUnits(config.InputUnits)
	if len(unknownUnits) > 0 {
		fmt.Printf("found input unit conflict: %v\n", unknownUnits)
		os.Exit(0)
	}
	config.OutputUnits = mergeEmpty(config.InputUnits, config.OutputUnits)
	config.OutputUnits, unknownUnits = checkUnits(config.OutputUnits)
	if len(unknownUnits) > 0 {
		fmt.Printf("found output unit conflict: %v\n", unknownUnits)
		os.Exit(0)
	}

	if len(config.ModelTable) > 0 {
		if len(config.Models) > 0 {
			fmt.Printf("Simultaneous table model listing and direct model specification not supported\n")
			os.Exit(0)
		}
		file, err := os.Open(config.ModelTable)
		if err != nil {
			fmt.Printf("Model table opening error: %v\n", err)
			os.Exit(0)
		}
		defer file.Close()

		models := []*ModelParameters{}
		if err := gocsv.UnmarshalFile(file, &models); err != nil {
			fmt.Printf("Model table reading error: %v\n", err)
			os.Exit(0)
		}

		modelFieldNames, aliases := utils.GetStructFieldNamesAndAliases(*models[0])
		header, err := utils.GetHeader(config.ModelTable)
		if err != nil {
			fmt.Printf("Model table error: %v\n", err)
			os.Exit(0)
		}
		for i := range header {
			name, _, _ := strings.Cut(header[i], "(")
			header[i] = strings.Trim(name, " ")
		}
		fieldsRead := []string{} //utils.FullIntersection(modelFieldNames, header)
		for n := range modelFieldNames {
			for h := range header {
				if header[h] == modelFieldNames[n] || (header[h] == aliases[n] && aliases[n] != "") {
					fieldsRead = append(fieldsRead, modelFieldNames[n])
				}
			}
		}

		groupname := utils.GetFilename(config.ModelTable)
		config.Models = make(map[string]ModelParameters, len(models))
		for line := range models {
			modelName := groupname + "_l" + strconv.Itoa(line+1)
			config.Models[modelName] = *models[line]
			for fieldName := range fieldsRead {
				// config.Models[modelName]._isFieldDefined[fieldsRead[fieldName]] = struct{}{}
				config.isDefinedMap[strings.Join([]string{"Models", modelName, fieldsRead[fieldName]}, "#")] = struct{}{}
			}
		}
	} else {
		if len(config.Models) == 0 {
			fmt.Println("No models provided")
			os.Exit(0)
		}
	}

	for modelName := range config.Models {
		model := config.Models[modelName]
		model._prototypeName = modelName
		config.Models[modelName] = model
	}

	if len(config.Options) > 0 {
		temporary := config.Models
		config.Models = make(map[string]ModelParameters, len(temporary))
		for modelName := range temporary {
			config.Models["/"+modelName] = temporary[modelName]
		}

		groups := make([]string, 0, len(config.Options))
		for group := range config.Options {
			groups = append(groups, group)
		}
		slices.Sort(groups)
		slices.Reverse(groups)
		for g := range groups {
			if len(config.Options[groups[g]]) == 0 {
				fmt.Printf("No options provided for parameter group: [%s]\n", groups[g])
				os.Exit(0)
			}
			draft := config.Models
			config.Models = make(map[string]ModelParameters, len(draft)*len(config.Options[groups[g]]))
			for draftName := range draft {
				for choice := range config.Options[groups[g]] {
					modelName := groups[g] + "-" + choice + "_" + draftName
					mp := draft[draftName]
					choiceElements := config.Options[groups[g]][choice]
					for ce := range choiceElements {
						choiceElementDescription := strings.Split(choiceElements[ce], "=")
						if len(choiceElementDescription) != 2 {
							fmt.Printf("Invalid optional parameter specification: [%s]\n", choiceElements[ce])
							os.Exit(0)
						}

						parameter := strings.Trim(choiceElementDescription[0], " ")
						valueStr := strings.Trim(choiceElementDescription[1], " ")

						modelField := reflect.ValueOf(&mp).Elem().FieldByName(parameter)
						if !modelField.IsValid() {
							fmt.Printf("Unknown parameter in optional specification: [%s]\n", parameter)
							os.Exit(0)
						}
						var valueReflect reflect.Value
						var err error
						switch modelField.Kind() {
						case reflect.Invalid:
							fmt.Printf("Unknown parameter in optional specification: %s\n", parameter)
							os.Exit(0)
						case reflect.Bool:
							var value bool
							value, err = strconv.ParseBool(valueStr)
							valueReflect = reflect.ValueOf(value)
						case reflect.Float64:
							var value float64
							value, err = strconv.ParseFloat(valueStr, 64)
							valueReflect = reflect.ValueOf(value)
						case reflect.Int:
							var value int
							value, err = strconv.Atoi(valueStr)
							valueReflect = reflect.ValueOf(value)
						case reflect.String:
							valueReflect = reflect.ValueOf(valueStr)
						default:
							fmt.Printf("Unsettable parameter in optional specification: %s\n", parameter)
							os.Exit(0)
						}
						if err != nil {
							fmt.Printf("Unable to parse value %s for field %s\n", valueStr, parameter)
							os.Exit(0)
						}
						modelField.Set(valueReflect)
						config.isDefinedMap[strings.Join([]string{"Models", modelName, parameter}, "#")] = struct{}{}
					}

					for def := range config.isDefinedMap {
						descr := strings.Split(def, "#")
						if len(descr) == 3 && descr[1] == draftName {
							descr[1] = modelName
							_, exist := config.isDefinedMap[strings.Join(descr, "#")]
							if exist {
								fmt.Printf("Redefinition of %s for %s, was in %s\n", descr[2], modelName, draftName)
								os.Exit(0)
							}
							config.isDefinedMap[strings.Join(descr, "#")] = struct{}{}
						}
					}
					config.Models[modelName] = mp
				}
			}
		}
	}

	if config.OutputDir != "" && config.OutputDir != "." {
		config.OutputDir = strings.TrimSuffix(config.OutputDir, "/") + "/" + utils.GetFilename(*flags.ConfigFileNamePointer) + "/"
		fmt.Printf("OUTPUT DIR: %s\n", config.OutputDir)
		os.MkdirAll(config.OutputDir, 0750)
	} else {
		panic(fmt.Errorf("output path not specified"))
	}

	undecoded := meta.Undecoded()
	if len(undecoded) > 0 {
		fmt.Println("Found unknown fields: [", undecoded, "]. Aborting.")
		os.Exit(0)
	}

	config._verbose = *flags.Verbose
	config._threads = *flags.Threads

	return config, meta
}

type ModelParameters struct {
	CrossSections                           string
	ElasticScatteringMode                   string
	InelasticScatteringMode                 string
	IonizationEnergySharing                 string
	IonizationScatteringMode                string
	Species                                 string
	ThresholdType                           string
	RequireCollisionRelativeMargin          map[string]float64
	GapLength                               float64 `csv:"L" units:"Length:1"`
	PressureGapLength                       float64 `csv:"pL" units:"Pressure:1,Length:1"`
	CathodeFallLength                       float64 `csv:"d" units:"Length:1"`
	PressureCathodeFallLength               float64 `csv:"pd" units:"Pressure:1,Length:1"`
	CathodeFallPotential                    float64 `csv:"Voltage"`
	CathodeCurrentDensity                   float64 `csv:"j" units:"Current:1,Length:-2"`
	CathodeCurrentDensityPerPressureSquared float64 `csv:"j/p2" units:"Current:1,Length:-2,Pressure:-2"`
	CathodeCurrent                          float64 `csv:"I" units:"Current:1"`
	SecondaryEmissionCoefficient            float64
	UseDonkosSEC                            bool
	UniformField                            bool
	ConstEField                             float64 `units:"Voltage:1,Length:-1"`
	Temperature                             float64
	Pressure                                float64 `csv:"p" units:"Pressure:1"`
	CathodeRadius                           float64 `units:"Length:1"`
	TubeRadius                              float64 `units:"Length:1"`

	SourceIntegralRelativeMargin float64

	// AmbipolarDiffusionCoefficient float64
	SlowElectronTemperature float64

	EnergyStep, AngleStep                          float64 // [eV]
	EnergyDiscretizationStep, MuDiscretizationStep float64
	NElectrons                                     int
	AddByNElectrons                                int
	MakeDir                                        bool
	ParallelPlaneHollowCathode                     bool
	Volumetric                                     bool
	AnodeBackscatteringCoefficient0                float64
	AnodeBackscatteringCoefficientB                float64
	AnodeBackscatteringEnergyLossFraction          float64

	TubeBackscatteringCoefficient0       float64
	TubeBackscatteringCoefficientB       float64
	TubeBackscatteringEnergyLossFraction float64

	CalculateSecondaryEmissionCoefficient bool
	CalculateCurrentDensity               bool
	CalculateVoltage                      bool
	NoDiffusionLoss                       bool
	SimplifiedDiffusionScale              bool
	CathodeFallLengthPrecision            float64

	CountNulls            bool
	CalculateDistribution bool
	DebugOutput           bool

	SupressSpinner bool

	GasDensity            float64
	_crossSections        *lxgata.Collisions
	_outputUnits          UnitConfig
	_verbose              bool
	_threads              int
	_prototypeName        string
	_lowerEnergyThreshold float64
	// _isFieldDefined map[string]struct{}
	_elasticScatteringMode   lxgata.ScatteringMode
	_inelasticScatteringMode lxgata.ScatteringMode
	_ionizationEnergySharing IonizationEnergySharing
	_ionizationScattering    IonizationScattering
}

func (p *ModelParameters) CrossSectionsData() *lxgata.Collisions {
	return p._crossSections
}

func (p *ModelParameters) SetCrossSectionsData(cd *lxgata.Collisions) {
	p._crossSections = cd
}

func (p *ModelParameters) OutputUnits() UnitConfig {
	return p._outputUnits
}

func (p *ModelParameters) OutputUnit(uc UnitClass) string {
	return p._outputUnits.GetForClass(uc)
}

func (p *ModelParameters) SetOutputUnits(uconf UnitConfig) {
	p._outputUnits = uconf
}

func (p *ModelParameters) Verbose() bool {
	return p._verbose
}

func (p ModelParameters) PrototypeName() string {
	return p._prototypeName
}

func (p *ModelParameters) Threads() int {
	return p._threads
}

func (p *ModelParameters) GetScatteringModes() (elastic, inelastic lxgata.ScatteringMode) {
	return p._elasticScatteringMode, p._inelasticScatteringMode
}

func (p *ModelParameters) GetIonizationEnergySharingMode() IonizationEnergySharing {
	return p._ionizationEnergySharing
}

func (p *ModelParameters) GetIonizationScatteringMode() IonizationScattering {
	return p._ionizationScattering
}

func (p *ModelParameters) SimulationLength() float64 {
	if p.ParallelPlaneHollowCathode {
		return p.GapLength / 2
	} else {
		return p.GapLength
	}
}

func (p *ModelParameters) LowerEnergyThreshold() float64 {
	if p._lowerEnergyThreshold == 0 {
		p._lowerEnergyThreshold = p.CrossSectionsData().MinThresholdOfKind(lxgata.CollisionType(strings.ToUpper(p.ThresholdType)))
	}
	return p._lowerEnergyThreshold
}

func ReducedFieldAtCathode(V, dc, gasDensity float64) float64 {
	return 2 * V / (dc * gasDensity * constants.Townsend)
}

func ReducedFieldMidSheath(V, dc, gasDensity float64) float64 {
	return V / (dc * gasDensity * constants.Townsend)
}

func (p *ModelParameters) ReducedFieldAtCathode() float64 {
	return ReducedFieldAtCathode(p.CathodeFallPotential, p.CathodeFallLength, p.GasDensity)
}
func (p *ModelParameters) ReducedFieldMidSheath() float64 {
	return ReducedFieldMidSheath(p.CathodeFallPotential, p.CathodeFallLength, p.GasDensity)
}

var defaultValues = map[string]any{ // in SI-eV
	"ThresholdType":                         "Ionization",
	"IonizationEnergySharing":               "UniformRandom",
	"IonizationScatteringMode":              "Boeuf",
	"Pressure":                              101325. / 760., //[Pa]
	"UniformField":                          false,
	"ConstEField":                           -100., //[V/m]
	"Temperature":                           300.,  //[K]
	"CathodeFallLengthPrecision":            1e-4,  //[m]
	"EnergyStep":                            0.01,  //[eV]
	"AngleStep":                             45.,
	"NElectrons":                            1000,
	"MakeDir":                               true,
	"ParallelPlaneHollowCathode":            false,
	"CalculateSecondaryEmissionCoefficient": false,
	"CalculateCurrentDensity":               false,
	"Volumetric":                            false,
	"CountNulls":                            false,
	"EnergyDiscretizationStep":              0.1,
	"MuDiscretizationStep":                  1 / 90.,
	"SourceIntegralRelativeMargin":          0.05,
}

var fieldsXor = map[string][]string{
	"CalculateSecondaryEmissionCoefficient": {"CathodeFallLength", "CalculateCurrentDensity", "CalculateVoltage"},
	"CalculateCurrentDensity":               {"CathodeFallLength", "CalculateSecondaryEmissionCoefficient", "CathodeCurrentDensity", "CathodeCurrent", "CalculateVoltage"},
	"CalculateVoltage":                      {"CathodeFallPotential", "CalculateSecondaryEmissionCoefficient", "CalculateCurrentDensity"},
	"CathodeFallLength":                     {"CalculateSecondaryEmissionCoefficient"},
	"PressureGapLength":                     {"Pressure"},
	"Pressure":                              {"PressureGapLength"},
}

var fieldsAnd = map[string][]string{
	"Volumetric":                            {"CathodeRadius"},
	"CathodeCurrent":                        {"CathodeRadius"},
	"CalculateSecondaryEmissionCoefficient": {"Species"},
	"CalculateCurrentDensity":               {"Species"},
	"CalculateVoltage":                      {"Species"},
	"PressureGapLength":                     {"GapLength"},
	"PressureCathodeFallLength":             {"Pressure"},
}
var fieldsDerivable map[string][]string = map[string][]string{
	"CathodeCurrent":                          {"CathodeCurrentDensity"},
	"PressureGapLength":                       {"Pressure"},
	"PressureCathodeFallLength":               {"CathodeFallLength"},
	"CathodeCurrentDensityPerPressureSquared": {"CathodeCurrentDensity"},
}

var calculableFields = map[string]func(
	*ModelParameters,
	[]string,
) []string{
	"CathodeCurrent": func(mp *ModelParameters, definedFields []string) []string {
		if slices.Contains(definedFields, "CathodeRadius") {
			area := math.Pi * mp.CathodeRadius * mp.CathodeRadius
			mp.CathodeCurrentDensity = mp.CathodeCurrent / area
			return []string{"CathodeCurrentDensity"}
		}
		// fmt.Printf("field 'CathodeRadius' not found: required by CathodeCurrentDensity calculation from CathodeCurrent\n")
		return nil
	},
	"CathodeCurrentDensityPerPressureSquared": func(mp *ModelParameters, definedFields []string) []string {
		if slices.Contains(definedFields, "Pressure") {
			mp.CathodeCurrentDensity = mp.CathodeCurrentDensityPerPressureSquared * (mp.Pressure * mp.Pressure)
			return []string{"CathodeCurrentDensity"}
		}
		// fmt.Printf("field 'Pressure' not found: required by CathodeCurrentDensity calculation from CathodeCurrent\n")
		return nil
	},
	"PressureGapLength": func(mp *ModelParameters, definedFields []string) []string {
		if slices.Contains(definedFields, "GapLength") {
			mp.Pressure = mp.PressureGapLength / mp.GapLength
			return []string{"Pressure"}
		}
		// fmt.Printf("field 'GapLength' not found: required by Pressure calculation from PressureGapLength\n")
		return nil
	},
	"PressureCathodeFallLength": func(mp *ModelParameters, definedFields []string) []string {
		if slices.Contains(definedFields, "Pressure") {
			mp.CathodeFallLength = mp.PressureCathodeFallLength / mp.Pressure
			return []string{"CathodeFallLength"}
		}
		return nil
	},
}

func AnyToSIeV(target any, units UnitConfig, direct bool) {
	targetReflect := reflect.ValueOf(target).Elem()
	targetType := targetReflect.Type()
	for i := 0; i < targetReflect.NumField(); i++ {
		field := targetType.Field(i)
		value := targetReflect.Field(i)
		if field.Anonymous && field.Type.Kind() == reflect.Struct {
			embeddedValue := value
			embeddedType := embeddedValue.Type()
			for j := 0; j < embeddedValue.NumField(); j++ {
				embeddedField := embeddedType.Field(j)
				if embeddedValue.FieldByName(embeddedField.Name).CanFloat() {
					fieldUnits := embeddedField.Tag.Get("units")
					unitElements := GetUnitElements(fieldUnits)
					if len(unitElements) != 0 {
						value := embeddedValue.FieldByName(embeddedField.Name).Float()
						value = FieldToSIeV(value, unitElements, units, direct)
						embeddedValue.FieldByName(embeddedField.Name).SetFloat(value)
					}
				}
			}
		} else {
			fieldUnits := field.Tag.Get("units")
			unitElements := GetUnitElements(fieldUnits)
			if len(unitElements) != 0 {
				if targetReflect.FieldByName(field.Name).CanFloat() {
					value := targetReflect.FieldByName(field.Name).Float()
					value = FieldToSIeV(value, unitElements, units, direct)
					targetReflect.FieldByName(field.Name).SetFloat(value)
				}
			}
		}
	}

	// for name := range parameterNames {
	// 	if targetReflect.FieldByName(parameterNames[name]).CanFloat() {
	// 		value := targetReflect.FieldByName(parameterNames[name]).Float()
	// 		value = FieldToSIeV(value, valueUnits[parameterNames[name]], units, direct)
	// 		targetReflect.FieldByName(parameterNames[name]).SetFloat(value)
	// 	}
	// }
}

func MakeHeader(target any, units UnitConfig) (header []string) {
	targetReflect := reflect.Indirect(reflect.ValueOf(target))
	targetType := targetReflect.Type()
	for i := 0; i < targetReflect.NumField(); i++ {
		field := targetType.Field(i)
		value := targetReflect.Field(i)
		if field.Anonymous && field.Type.Kind() == reflect.Struct {
			embeddedValue := value
			embeddedType := embeddedValue.Type()
			for j := 0; j < embeddedValue.NumField(); j++ {
				embeddedField := embeddedType.Field(j)
				tag := embeddedField.Tag.Get("csv")
				if tag == "-" {
					break
				}
				var csvFieldName string
				if tag != "" {
					csvFieldName = tag
				} else {
					csvFieldName = embeddedField.Name
				}
				unitElements := GetUnitElements(embeddedField.Tag.Get("units"))
				if len(unitElements) != 0 {
					header = append(header, csvFieldName+" "+unitElements.String(units))
				} else {
					header = append(header, csvFieldName)
				}
			}
		} else {
			tag := field.Tag.Get("csv")
			if tag == "-" {
				break
			}
			var csvFieldName string
			if tag != "" {
				csvFieldName = tag
			} else {
				csvFieldName = field.Name
			}
			fieldUnits := field.Tag.Get("units")
			unitElements := GetUnitElements(fieldUnits)
			if len(unitElements) != 0 {
				header = append(header, csvFieldName+" "+unitElements.String(units))
			} else {
				header = append(header, csvFieldName)
			}
		}
	}

	return header
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

	AnyToSIeV(modelConfig, config.InputUnits, true)

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

	if !slices.Contains(discoveredParameters, "TubeRadius") {
		modelConfig.TubeRadius = modelConfig.CathodeRadius
	}

	scatteringMode := map[string]lxgata.ScatteringMode{
		"Born":            lxgata.Born,
		"ScreenedCoulomb": lxgata.Coulomb,
		"Isotropic":       lxgata.Isotropic,
	}
	if elsm, ok := scatteringMode[modelConfig.ElasticScatteringMode]; !ok {
		fmt.Printf("Wrong elastic scattering mode: %v\n", modelConfig.ElasticScatteringMode)
		allGood = false
	} else {
		modelConfig._elasticScatteringMode = elsm
	}
	if insm, ok := scatteringMode[modelConfig.InelasticScatteringMode]; !ok {
		fmt.Printf("Wrong inelastic scattering mode: %v\n", modelConfig.InelasticScatteringMode)
		allGood = false
	} else {
		modelConfig._inelasticScatteringMode = insm
	}
	ionizationScatteringMode := map[string]IonizationScattering{
		"Boeuf":           Boeuf,
		"Isotropic":       Isotropic,
		"ScreenedCoulomb": Coulomb,
		"CoulombM":        CoulombM,
		"Born":            Born,
		"BornFL":          BornFullLoss,
		"Mixed":           Mixed,
		"MixedFL":         MixedFullLoss,
	}
	if izsm, ok := ionizationScatteringMode[modelConfig.IonizationScatteringMode]; !ok {
		fmt.Printf("Wrong ionization scattering mode: %v\n", modelConfig.IonizationScatteringMode)
		allGood = false
	} else {
		modelConfig._ionizationScattering = izsm
	}

	energySharingMode := map[string]IonizationEnergySharing{
		"Zero":          Zero,
		"UniformRandom": UniformRandom,
		"Equal":         Equal,
		"Opal":          Opal,
	}
	if esm, ok := energySharingMode[modelConfig.IonizationEnergySharing]; !ok {
		fmt.Printf("Wrong energy sharing mode: %v\n", modelConfig.IonizationEnergySharing)
		allGood = false
	} else {
		modelConfig._ionizationEnergySharing = esm
	}

	modelConfig._threads = config._threads
	modelConfig._verbose = config._verbose

	return allGood
}
