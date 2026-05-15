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
	"github.com/wildstyl3r/stmc/internal/messages"
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

const (
	Collisionless float64 = 1. / 3.
	ConstMobility float64 = 1. / 2.
	StrongField   float64 = 2. / 3.
	Linear        float64 = 1
)

type UnitConfig struct {
	L string `form:"label=Length,options=mm|cm|m" toml:"L,omitzero,omitempty"`
	T string `form:"hide" toml:"T,omitzero,omitempty"`
	P string `form:"label=Pressure,options=Torr|Pa|mbar|bar" toml:"P,omitzero,omitempty"`
	I string `form:"label=Current,options=mkA|mA|A" toml:"I,omitzero,omitempty"`
	E string `form:"label=Energy,options=eV|J" toml:"E,omitzero,omitempty"`
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
	OutputDir    string     `form:"label=Output path,widget=folder" toml:"OutputDir,omitzero,omitempty"`
	InputUnits   UnitConfig `form:"widget=row" toml:"InputUnits,omitzero,omitempty"`
	OutputUnits  UnitConfig `form:"widget=row" toml:"OutputUnits,omitzero,omitempty"`
	AddPotential float64    `form:"label=Potential offset" toml:"AddPotential,omitzero,omitempty"`

	ModelParameters `form:"label=Global defaults,sparse"`
	ModelTable      string                                 `form:"label=Models from file,widget=file" toml:"ModelTable,omitzero,omitempty"`
	Models          map[string]*ModelParameters            `form:"element=Model" toml:"Models,omitzero,omitempty"`
	Options         map[string]map[string]*ModelParameters `form:"element=Option Group,widget=options" toml:"Options,omitzero,omitempty"` //[]string
	isDefinedMap    map[string]struct{}
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

func LoadConfig(flags Flags, logger messages.Logger) (Config, toml.MetaData) {
	var config Config
	config.isDefinedMap = map[string]struct{}{}
	meta, err := toml.DecodeFile(strings.TrimSuffix(*flags.ConfigFileNamePointer, ".toml")+".toml", &config)
	if err != nil {
		// fmt.Println(err)
		// fmt.Fprintln(os.Stdout, err)
		panic(err)
	}

	var unknownUnits []string
	config.InputUnits, unknownUnits = checkUnits(config.InputUnits)
	if len(unknownUnits) > 0 {
		logger.Failure("found input unit conflict: %v\n", unknownUnits)
	}
	config.OutputUnits = mergeEmpty(config.InputUnits, config.OutputUnits)
	config.OutputUnits, unknownUnits = checkUnits(config.OutputUnits)
	if len(unknownUnits) > 0 {
		logger.Failure("found input unit conflict: %v\n", unknownUnits)
	}

	if len(config.ModelTable) > 0 {
		if len(config.Models) > 0 {
			logger.Failure("Simultaneous table model listing and direct model specification not supported\n")
		}
		file, err := os.Open(config.ModelTable)
		if err != nil {
			logger.Failure("Model table opening error: %v\n", err)
		}
		defer file.Close()

		models := []*ModelParameters{}
		if err := gocsv.UnmarshalFile(file, &models); err != nil {
			logger.Failure("Model table reading error: %v\n", err)
		}

		modelFieldNames, aliases := utils.GetStructFieldNamesAndAliases(*models[0])
		header, err := utils.GetHeader(config.ModelTable)
		if err != nil {
			logger.Failure("Model table error: %v\n", err)
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
		config.Models = make(map[string]*ModelParameters, len(models))
		for line := range models {
			modelName := groupname + "_l" + strconv.Itoa(line+1)
			baseModel := *models[line]
			config.Models[modelName] = &baseModel
			for fieldName := range fieldsRead {
				// config.Models[modelName]._isFieldDefined[fieldsRead[fieldName]] = struct{}{}
				config.isDefinedMap[strings.Join([]string{"Models", modelName, fieldsRead[fieldName]}, "#")] = struct{}{}
			}
		}
	} else {
		if len(config.Models) == 0 {
			logger.Failure("No models provided")
		}
	}

	for modelName := range config.Models {
		config.Models[modelName]._prototypeName = modelName
		// model._prototypeName = modelName
		// config.Models[modelName] = model
	}

	if len(config.Options) > 0 {
		temporary := config.Models
		config.Models = make(map[string]*ModelParameters, len(temporary))
		for modelName := range temporary {
			t := *temporary[modelName]
			config.Models["/"+modelName] = &t
		}

		groups := make([]string, 0, len(config.Options))
		for group := range config.Options {
			groups = append(groups, group)
		}
		slices.Sort(groups)
		slices.Reverse(groups)
		for g := range groups {
			if len(config.Options[groups[g]]) == 0 {
				logger.Failure("No options provided for parameter group: [%s]\n", groups[g])
			}
			draft := config.Models
			config.Models = make(map[string]*ModelParameters, len(draft)*len(config.Options[groups[g]]))
			for draftName := range draft {
				for choice := range config.Options[groups[g]] {
					modelName := groups[g] + "-" + choice + "_" + draftName
					mp := *draft[draftName]
					choiceParameters := config.Options[groups[g]][choice]
					choiceReflect := reflect.TypeOf(*choiceParameters)
					for i := 0; i < choiceReflect.NumField(); i++ {
						parameter := choiceReflect.Field(i).Name
						if meta.IsDefined("Options", groups[g], choice, parameter) {
							modelField := reflect.ValueOf(&mp).Elem().FieldByName(parameter)
							if !modelField.IsValid() {
								logger.Failure("Unknown parameter in optional specification: [%s]\n", parameter)
							}
							dummyField := reflect.ValueOf(choiceParameters).Elem().FieldByName(parameter)
							modelField.Set(dummyField)
							config.isDefinedMap[strings.Join([]string{"Models", modelName, parameter}, "#")] = struct{}{}
						}
					}
					// for ce := range choiceParameters {
					// 	dummy := ModelParameters{}
					// 	_, dummyError := toml.Decode(choiceParameters[ce], &dummy)
					// 	if dummyError != nil {
					// 		fmt.Printf("error decoding option: %v", dummyError)
					// 		os.Exit(0)
					// 	}
					// 	choiceElementDescription := strings.Split(choiceParameters[ce], "=")

					// 	parameter := strings.Trim(choiceElementDescription[0], " ")
					// 	// valueStr := strings.Trim(choiceElementDescription[1], " ")

					// }

					for def := range config.isDefinedMap {
						descr := strings.Split(def, "#")
						if len(descr) == 3 && descr[1] == draftName {
							descr[1] = modelName
							_, exist := config.isDefinedMap[strings.Join(descr, "#")]
							if exist {
								logger.Failure("Redefinition of %s for %s, was in %s\n", descr[2], modelName, draftName)
							}
							config.isDefinedMap[strings.Join(descr, "#")] = struct{}{}
						}
					}
					config.Models[modelName] = &mp
				}
			}
		}
	}

	if config.OutputDir != "" && config.OutputDir != "." {
		config.OutputDir = strings.TrimSuffix(config.OutputDir, "/") + "/" + utils.GetFilename(*flags.ConfigFileNamePointer) + "/"
		fmt.Printf("OUTPUT DIR: %s\n", config.OutputDir)
		os.MkdirAll(config.OutputDir, 0750)
	} else {
		logger.Failure("output path not specified")
	}

	undecoded := meta.Undecoded()
	if len(undecoded) > 0 {
		logger.Failure("Found unknown fields: [", undecoded, "]. Aborting.")
	}

	config._verbose = *flags.Verbose
	config._threads = *flags.Threads

	return config, meta
}

// type ChemicalComposition struct {
// 	ComponentsShare    map[string]float64
// 	DissociativeFactor int
// 	InitialShare float64
// }

type CalculationMode int

const (
	Unspecified CalculationMode = iota
	BasicSheathCalculation
	TownsendAlpha
	GammaCalculation
	VoltageCalculation
	CurrentCalculation
	IonMobilityCalculation
)

type EmissionMode int

const (
	Cosine EmissionMode = iota
	ForwardIsotropic
	Forward
)

type ModelParameters struct {
	SheathFieldModel                        string             `form:"label=Sheath field model,options=Linear|ConstMobility|Collisionless|StrongField" toml:"SheathFieldModel,omitzero,omitempty"`
	CrossSections                           string             `form:"label=Cross section file,widget=file" toml:"CrossSections,omitzero,omitempty"`
	CrossSectionsDStep                      float64            `form:"label=Cross section discretization step" toml:"CrossSectionsDStep,omitzero,omitempty"`
	ElasticScatteringMode                   string             `form:"label=Elastic scattering model,options=BornBethe|ScreenedCoulomb|Isotropic" toml:"ElasticScatteringMode,omitzero,omitempty"`
	InelasticScatteringMode                 string             `form:"label=Inelastic scattering model,options=BornBethe|ScreenedCoulomb|Isotropic" toml:"InelasticScatteringMode,omitzero,omitempty"`
	IonizationEnergySharing                 string             `form:"label=Ionization energy sharing,options=Opal|UniformRandom|Equal" toml:"IonizationEnergySharing,omitzero,omitempty"`
	IonizationScatteringMode                string             `form:"label=Ionization scattering model,options=Boeuf|Isotropic|ScreenedCoulomb|Born" toml:"IonizationScatteringMode,omitzero,omitempty"`
	Species                                 map[string]float64 `form:"label=Species,widget=floatmap,element=Species" toml:"Species,omitzero,omitempty"`
	ThresholdType                           string             `form:"label=Lower energy cutoff by threshold of type" toml:"ThresholdType,omitzero,omitempty"`
	LowerThresholdValue                     float64            `form:"label=Cutoff energy" toml:"LowerThresholdValue,omitzero,omitempty"`
	NoCutoff                                bool               `form:"label=Do not use lower energy cutoff" toml:"NoCutoff,omitzero,omitempty"`
	RequireCollisionRelativeMargin          map[string]float64 `form:"label=Precision for collisions of interest,widget=floatmap,element=Process Wildcard" toml:"RequireCollisionRelativeMargin,omitzero,omitempty"`
	GapLength                               float64            `form:"label=Gap length" csv:"L" units:"Length:1" toml:"GapLength,omitzero,omitempty"`
	PressureGapLength                       float64            `form:"label=Barometric gap length" csv:"pL" units:"Pressure:1,Length:1" toml:"PressureGapLength,omitzero,omitempty"`
	CathodeFallLength                       float64            `form:"label=Sheath thickness" csv:"d" units:"Length:1" toml:"CathodeFallLength,omitzero,omitempty"`
	PressureCathodeFallLength               float64            `form:"label=Barometric sheath thickness" csv:"pd" units:"Pressure:1,Length:1" toml:"PressureCathodeFallLength,omitzero,omitempty"`
	CathodeFallPotential                    float64            `form:"label=Cathode fall voltage" csv:"Voltage" toml:"CathodeFallPotential,omitzero,omitempty"`
	CathodeCurrentDensity                   float64            `form:"label=Cathode j" csv:"j" units:"Current:1,Length:-2" toml:"CathodeCurrentDensity,omitzero,omitempty"`
	CathodeCurrentDensityPerPressureSquared float64            `form:"label=Cathode j/p^2" csv:"j/p2" units:"Current:1,Length:-2,Pressure:-2" toml:"CathodeCurrentDensityPerPressureSquared,omitzero,omitempty"`
	CathodeCurrentPerPressureSquared        float64            `form:"label=Cathode I/p^2" csv:"I/p2" units:"Current:1,Pressure:-2" toml:"CathodeCurrentPerPressureSquared,omitzero,omitempty"`
	CathodeCurrent                          float64            `form:"label=Cathode I" csv:"I" units:"Current:1" toml:"CathodeCurrent,omitzero,omitempty"`
	SecondaryEmissionCoefficient            float64            `form:"label=Prescribed gamma" toml:"SecondaryEmissionCoefficient,omitzero,omitempty"`
	UseDonkosSEC                            bool               `toml:"UseDonkosSEC,omitzero,omitempty"`
	UniformField                            bool               `toml:"UniformField,omitzero,omitempty"`
	RandomizedWallLoss                      bool               `toml:"RandomizedWallLoss,omitzero,omitempty"`
	ConstEField                             float64            `form:"label=Const electric field in nglow" units:"Voltage:1,Length:-1" toml:"ConstEField,omitzero,omitempty"`
	Temperature                             float64            `toml:"Temperature,omitzero,omitempty"`
	Pressure                                float64            `toml:"Pressure,omitzero,omitempty" csv:"p" units:"Pressure:1"`
	CathodeRadius                           float64            `toml:"CathodeRadius,omitzero,omitempty" units:"Length:1"`
	TubeRadius                              float64            `toml:"TubeRadius,omitzero,omitempty" units:"Length:1"`

	SourceIntegralRelativeMargin float64 `form:"label=Relative margin for gamma and int[S(x)]" toml:"SourceIntegralRelativeMargin,omitzero,omitempty"`

	// AmbipolarDiffusionCoefficient float64
	SlowElectronTemperature float64 `form:"label=Slow electron temperature" toml:"SlowElectronTemperature,omitzero,omitempty"`

	EnergyStep                            float64 `form:"label=Energy step" toml:"EnergyStep,omitzero,omitempty"` // [eV]
	MuDiscretizationStep                  float64 `form:"hide" toml:"MuDiscretizationStep,omitzero,omitempty"`
	NParticles                            int     `form:"label=Starting electron number" toml:"NParticles,omitzero,omitempty"`
	AddByNElectrons                       int     `form:"label=Additional electron number" toml:"AddByNElectrons,omitzero,omitempty"`
	ParallelPlaneHollowCathode            bool    `toml:"ParallelPlaneHollowCathode,omitzero,omitempty"`
	Volumetric                            bool    `toml:"Volumetric,omitzero,omitempty"`
	AnodeBackscatteringCoefficient0       float64 `toml:"AnodeBackscatteringCoefficient0,omitzero,omitempty"`
	AnodeBackscatteringCoefficientB       float64 `toml:"AnodeBackscatteringCoefficientB,omitzero,omitempty"`
	AnodeBackscatteringEnergyLossFraction float64 `toml:"AnodeBackscatteringEnergyLossFraction,omitzero,omitempty"`

	TubeBackscatteringCoefficient0       float64 `toml:"TubeBackscatteringCoefficient0,omitzero,omitempty"`
	TubeBackscatteringCoefficientB       float64 `toml:"TubeBackscatteringCoefficientB,omitzero,omitempty"`
	TubeBackscatteringEnergyLossFraction float64 `toml:"TubeBackscatteringEnergyLossFraction,omitzero,omitempty"`

	ReturningElectrons bool `form:"label=Cutoff by cathode return-ability" toml:"ReturningElectrons,omitzero,omitempty"`

	// CalculateSecondaryEmissionCoefficient bool
	// CalculateCurrentDensity               bool
	// CalculateVoltage           bool
	CalculationMode            string  `form:"label=Calculation mode,options=Basic|Current|Gamma|Voltage" toml:"CalculationMode,omitzero,omitempty"` //CalculationMode
	NoDiffusionLoss            bool    `toml:"NoDiffusionLoss,omitzero,omitempty"`
	SimplifiedDiffusionScale   bool    `toml:"SimplifiedDiffusionScale,omitzero,omitempty"`
	CathodeFallLengthPrecision float64 `toml:"CathodeFallLengthPrecision,omitzero,omitempty"`

	CountCollisions []string `toml:"StoreCollisions,omitzero,omitempty"`

	EmissionType string `form:"label=Emission type,options=Cosine,ForwardIsotropic,Forward" toml:"EmissionType,omitzero,omitempty"`

	CalculateDistribution   bool `toml:"CalculateDistribution,omitzero,omitempty"`
	DebugOutput             bool `toml:"DebugOutput,omitzero,omitempty"`
	TrajectoryAveragedRates bool `form:"hide" toml:"TrajectoryAveragedRates,omitzero,omitempty"`
	EnergyDeposition        bool `toml:"EnergyDeposition,omitzero,omitempty"`
	MeanFreePath            bool `toml:"MeanFreePath,omitzero,omitempty"`
	DarkDischarge           bool `toml:"DarkDischarge,omitzero,omitempty"`

	SupressSpinner bool `form:"hide" toml:"-"`

	GasDensity            float64 `form:"hide" toml:"-"`
	_crossSections        *lxgata.Collisions
	_outputUnits          UnitConfig
	_verbose              bool
	_superelastics        bool
	_threads              int
	_prototypeName        string
	_lowerEnergyThreshold float64
	_mixtureParameters    map[string]lxgata.Species
	// _isFieldDefined map[string]struct{}
	_elasticScatteringMode   lxgata.ScatteringMode
	_inelasticScatteringMode lxgata.ScatteringMode
	_ionizationEnergySharing IonizationEnergySharing
	_ionizationScattering    IonizationScattering
	_calculationMode         CalculationMode
	_emissionMode            EmissionMode
	_sheathFieldModel        float64
	_countCollisions         map[lxgata.CollisionType]struct{}
}

func (p *ModelParameters) CrossSectionsData() *lxgata.Collisions {
	return p._crossSections
}

func (p *ModelParameters) SetCrossSectionsData(cd *lxgata.Collisions) {
	p._crossSections = cd
	for process := range cd.Processes {
		if cd.Processes[process].Type == lxgata.DEEXCITATION {
			p._superelastics = true
			return
		}
	}
}

func (p *ModelParameters) GetCalculationMode() CalculationMode {
	return p._calculationMode
}

func (p *ModelParameters) GetEmissionMode() EmissionMode {
	return p._emissionMode
}

func (p *ModelParameters) HasSuperelastics() bool {
	return p._superelastics
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

func (p *ModelParameters) GetSheathFieldPower() float64 {
	return p._sheathFieldModel
}

func (p *ModelParameters) GetMixtureParameters() map[string]lxgata.Species {
	return p._mixtureParameters
}

func (p *ModelParameters) GetSpeciesString() string {
	elements := []string{}
	for key := range p.Species {
		elements = append(elements, key)
	}
	slices.Sort(elements)
	return strings.Join(elements, " ")
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

func (p *ModelParameters) GetCollisionTypesToStore() map[lxgata.CollisionType]struct{} {
	return p._countCollisions
}

func (p *ModelParameters) SimulationLength() float64 {
	if p.ParallelPlaneHollowCathode {
		return p.GapLength / 2
	} else {
		return p.GapLength
	}
}

func (p *ModelParameters) LowerEnergyThreshold() float64 {
	if p.NoCutoff {
		return 0
	}
	if p._lowerEnergyThreshold == 0 {
		if p.LowerThresholdValue != 0 {
			p._lowerEnergyThreshold = p.LowerThresholdValue
			return p._lowerEnergyThreshold
		}
		p._lowerEnergyThreshold = p.CrossSectionsData().MinThresholdOfKind(lxgata.CollisionType(strings.ToUpper(p.ThresholdType))) - 0.5
		for process := range p.CrossSectionsData().Processes {
			if p.CrossSectionsData().Processes[process].Type == lxgata.DEEXCITATION {
				p._lowerEnergyThreshold = -1
				return -1
			}
		}
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
	"CountCollisions":            []string{"IONIZATION"},
	"CrossSectionsDStep":         0.001,
	"SheathFieldModel":           "Linear",
	"ThresholdType":              "Ionization",
	"IonizationEnergySharing":    "UniformRandom",
	"EmissionType":               "Cosine",
	"IonizationScatteringMode":   "Boeuf",
	"Pressure":                   101325. / 760., //[Pa]
	"UniformField":               false,
	"ConstEField":                -100., //[V/m]
	"Temperature":                300.,  //[K]
	"CathodeFallLengthPrecision": 1e-4,  //[m]
	"EnergyStep":                 0.01,  //[eV]
	"NParticles":                 1000,
	"ParallelPlaneHollowCathode": false,
	"CalculationMode":            "Basic",
	// "CalculateSecondaryEmissionCoefficient": false,
	// "CalculateCurrentDensity":               false,
	"Volumetric":                   false,
	"MuDiscretizationStep":         1 / 50.,
	"SourceIntegralRelativeMargin": 0.005,
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
	"CathodeCurrentPerPressureSquared":        {"CathodeCurrentDensity"},
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
	"CathodeCurrentPerPressureSquared": func(mp *ModelParameters, definedFields []string) []string {
		if slices.Contains(definedFields, "Pressure") {
			mp.CathodeCurrentDensity = mp.CathodeCurrentPerPressureSquared * (mp.Pressure * mp.Pressure)
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

func embeddedNames(target reflect.Value, units UnitConfig) (names []string) {
	targetReflect := reflect.Indirect(target)
	targetType := targetReflect.Type()
	for i := range targetReflect.NumField() {
		field := targetType.Field(i)
		value := targetReflect.Field(i)
		if field.Anonymous && field.Type.Kind() == reflect.Struct {
			names = append(names, embeddedNames(value, units)...)
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
				names = append(names, csvFieldName+" "+unitElements.String(units))
			} else {
				names = append(names, csvFieldName)
			}
		}
	}
	return names
}

func MakeHeader(target any, units UnitConfig) (header []string) {
	return embeddedNames(reflect.ValueOf(target), units)
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

func (modelConfig *ModelParameters) CheckAndUnify(modelName string, config *Config, meta *toml.MetaData, logger messages.Logger) bool {

	globalAmbiguities, globalMissingDeps := config.checkFieldProblems([]string{}, meta, config)
	localAmbiguities, localMissingDeps := modelConfig.checkFieldProblems([]string{"Models", modelName}, meta, config)
	if len(globalAmbiguities) > 0 {
		logger.Info("unable to load config: found global ambiguities \n%v\n", globalAmbiguities)
		return false
	}
	if len(localAmbiguities) > 0 {
		logger.Info("unable to load config: found model ambiguities \n%v\n", localAmbiguities)
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
		logger.Info("unable to load config: required dependent fields not found \n%v\n", missingIntersection)
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

	it := 0
	for fieldName := range defaultValues {
		it++
		if _, x := excludeFromLoadingDefaultOrOuter[fieldName]; !x && !slices.Contains(discoveredParameters, fieldName) {
			fv := modelConfigReflect.Elem().FieldByName(fieldName)
			if fv.IsValid() {
				if fv.CanSet() {
					fv.Set(reflect.ValueOf(defaultValues[fieldName]))
				} else {
					logger.Failure("cannot set %v\n", fv)
				}
			} else {
				logger.Failure("fv %s of %#v is invalid: %v, again: %v\n", fieldName, modelConfigReflect.Elem(), fv, modelConfigReflect.Elem().FieldByName(fieldName))
			}

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
	allGood := true
	status := ""

	for initialFieldName := range calculableFields {
		if slices.Contains(enabledParameters, initialFieldName) {
			allGood = false
			status += fmt.Sprintf("field conflict: calculable field %v in enabled fields %v\n", initialFieldName, enabledParameters)
		}
	}

	for i := range enabledParameters {
		for requirement := range fieldsAnd[enabledParameters[i]] {
			if !slices.Contains(enabledParameters, fieldsAnd[enabledParameters[i]][requirement]) {
				status += fmt.Sprintf("for parameter %s requirement %s not found\n", enabledParameters[i], fieldsAnd[enabledParameters[i]][requirement])
				allGood = false
			}
		}
		for conflict := range fieldsXor[enabledParameters[i]] {
			if slices.Contains(enabledParameters, fieldsXor[enabledParameters[i]][conflict]) {
				status += fmt.Sprintf("for parameter %s found conflicting parameter: %s\n", enabledParameters[i], fieldsXor[enabledParameters[i]][conflict])
				allGood = false
			}
		}
	}

	modelConfig.CathodeFallPotential += config.AddPotential
	modelConfig.GasDensity = modelConfig.Pressure / (constants.KBolzmann * modelConfig.Temperature)
	var conflict []string
	units, conflict := checkUnits(config.OutputUnits)
	if len(conflict) > 0 {
		status += fmt.Sprintf("found output unit conflict: %v\n Data will be saved in input units\n", conflict)
		modelConfig._outputUnits = config.InputUnits
	} else {
		modelConfig._outputUnits = units
	}

	if !slices.Contains(discoveredParameters, "TubeRadius") {
		modelConfig.TubeRadius = modelConfig.CathodeRadius
	}

	scatteringMode := map[string]lxgata.ScatteringMode{
		"BornBethe":       lxgata.Born,
		"ScreenedCoulomb": lxgata.Coulomb,
		"Isotropic":       lxgata.Isotropic,
	}
	if elsm, ok := scatteringMode[modelConfig.ElasticScatteringMode]; !ok {
		status += fmt.Sprintf("Wrong elastic scattering mode: %v\n", modelConfig.ElasticScatteringMode)
		allGood = false
	} else {
		modelConfig._elasticScatteringMode = elsm
	}
	if insm, ok := scatteringMode[modelConfig.InelasticScatteringMode]; !ok {
		status += fmt.Sprintf("Wrong inelastic scattering mode: %v\n", modelConfig.InelasticScatteringMode)
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
		status += fmt.Sprintf("Wrong ionization scattering mode: %v\n", modelConfig.IonizationScatteringMode)
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
		status += fmt.Sprintf("Wrong energy sharing mode: %v\n", modelConfig.IonizationEnergySharing)
		allGood = false
	} else {
		modelConfig._ionizationEnergySharing = esm
	}

	modelConfig._threads = config._threads
	modelConfig._verbose = config._verbose

	mixtureTotalWeight := 0.
	for _, weight := range modelConfig.Species {
		mixtureTotalWeight += weight
	}
	modelConfig._mixtureParameters = make(map[string]lxgata.Species)
	for species, weight := range modelConfig.Species {
		modelConfig._mixtureParameters[species] = lxgata.Species{
			ShareOfUnity: weight / mixtureTotalWeight,
			UParameter:   lxgata.Hartree,
		}
	}
	if mixtureTotalWeight == 0 {
		allGood = false
		status += "Gas mixture malformed: no species provided with non-zero weights\n"
	}

	calculationMode := map[string]CalculationMode{
		"Current": CurrentCalculation,
		"Gamma":   GammaCalculation,
		"Voltage": VoltageCalculation,
		"Basic":   BasicSheathCalculation,
		"Alpha":   TownsendAlpha,
	}
	modelConfig._calculationMode = calculationMode[modelConfig.CalculationMode]

	emissionMode := map[string]EmissionMode{
		"Cosine":           Cosine,
		"Forward":          Forward,
		"ForwardIsotropic": ForwardIsotropic,
	}
	modelConfig._emissionMode = emissionMode[modelConfig.EmissionType]

	sheathFieldModel := map[string]float64{
		"Linear":        Linear,
		"Collisionless": Collisionless,
		"ConstMobility": ConstMobility,
		"StrongField":   StrongField,
	}
	modelConfig._sheathFieldModel = sheathFieldModel[modelConfig.SheathFieldModel]

	modelConfig._countCollisions = make(map[lxgata.CollisionType]struct{})
	for _, t := range modelConfig.CountCollisions {
		modelConfig._countCollisions[lxgata.CollisionType(t)] = struct{}{}
	}

	if status != "" {
		logger.Info(status)
	}

	if modelConfig.GetCalculationMode() == TownsendAlpha {
		// modelConfig.NoCutoff = true
		modelConfig.DarkDischarge = true
	}

	return allGood
}

func StringifyMixture(mixture map[string]lxgata.Species) (s string) {
	s = "{"
	components := make([]struct {
		key string
		val float64
	}, 0, len(mixture))
	for key, val := range mixture {
		components = append(components, struct {
			key string
			val float64
		}{key, val.ShareOfUnity})

	}
	slices.SortFunc(components, func(a, b struct {
		key string
		val float64
	}) int {
		return strings.Compare(a.key, b.key)
	})
	for _, component := range components {
		s = s + component.key + ": " + strconv.FormatFloat(component.val, 'e', 6, 64) + "_"
	}
	s = s[:len(s)-1] + "}"
	return s
}
