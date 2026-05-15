package model

import (
	"fmt"
	"math"
	"math/rand/v2"
	"strings"
	"sync"

	"github.com/wildstyl3r/lxgata"
	"github.com/wildstyl3r/stmc/internal/config"
	"github.com/wildstyl3r/stmc/internal/constants"
	"github.com/wildstyl3r/stmc/internal/messages"
	"github.com/wildstyl3r/stmc/internal/utils"
)

type Model struct {
	Parameters                config.ModelParameters
	Vs                        float64 // cathode fall potential, V_sheath
	FieldSwitchPosition       float64
	FieldSwitchPotential      float64 // < 0
	EFieldScale               float64
	SheathPotentialCorrection float64
	//motionK, motionW          float64
	inverseCathodeFallLength float64
	inverseConstEField       float64
	inverseCathodePotential  float64
	sheathTimeOuterConstant  float64
	sheathTimeInnerConstant  float64
	constFieldTimeConstant   float64
	constFieldAcceleration   float64

	cathodeRadius2, tubeRadius2 float64

	// Va float64 // additional voltage to avoid numeric negative energy beyond cathode fall

	NumCells, NumCellsMu, NumCellsE int

	XStep float64

	lookUpPotentialAtCellNode   []float64
	LookupEnergy                []float64
	lookupMNumerator            []float64
	lookupCellNodeFromPotential []int

	DistributionXEMu                     [][][]float64
	IonizationCSmulGasDensityLookupECell []float64
	TrajAveragedIonization               []utils.KahanSummable
	TimeAtCell                           []utils.AggregationF
	MeanFreePath                         []utils.AggregationF
	GlobalMeanFreePath                   utils.AggregationF
	CollisionAtCell                      map[lxgata.CollisionType][]utils.Aggregation
	DetailedCollisionAtCell              map[string][]utils.Aggregation
	MeanEnergy                           []utils.AggregationF
	IonizationsSumUpToCell               []utils.Aggregation

	IonizationCounters     []int
	AttachmentCounters     []int
	AnodeElectronCounter   int
	CathodeElectronCounter int

	TotalElectronsEmittedOnCathode int
	ElectronsReturned              utils.Aggregation
	IonizingCathodeElectrons       utils.Aggregation

	MeanElectronEnergyAtAnode float64

	DataHub *DataHubType
}

func (m *Model) CoreResult() utils.Result {
	electronsReturned, electronsReturnedMargin := m.ElectronsReturned.MeanWithErrorMargin(0.95)
	ionizingElectrons, ionizingElectronsMargin := m.IonizingCathodeElectrons.MeanWithErrorMargin(0.95)
	return utils.Result{
		ModelName:                      m.Parameters.PrototypeName(),
		ReducedFieldAtCathode:          m.ReducedFieldAtCathode(),
		GlobalMeanFreePath:             m.GlobalMeanFreePath.Mean(),
		Voltage:                        m.Parameters.CathodeFallPotential,
		Pressure:                       m.Parameters.Pressure,
		ElectronsReturned:              electronsReturned,
		ElectronsReturnedMargin:        electronsReturnedMargin,
		IonizingCathodeElectrons:       ionizingElectrons,
		IonizingCathodeElectronsMargin: ionizingElectronsMargin,
	}
}

func (m *Model) SheathResult() utils.SheathResult {
	return utils.SheathResult{
		Result:                     m.CoreResult(),
		ReducedFieldAtSheathCenter: m.ReducedFieldMidSheath(),
		CathodeFallLength:          m.Parameters.CathodeFallLength,
		PressureCathodeFallLength:  m.Parameters.CathodeFallLength * m.Parameters.Pressure,
	}
}

func (m *Model) InitVars() {
	m.Vs = math.Abs(m.Parameters.CathodeFallPotential)
	// m.Va = math.Abs(m.Parameters.ConstEField * (m.Parameters.SimulationLength() - m.Parameters.CathodeFallLength))
	m.inverseCathodeFallLength = 1. / m.Parameters.CathodeFallLength
	m.inverseConstEField = 1. / m.Parameters.ConstEField
	m.inverseCathodePotential = 1. / m.Parameters.CathodeFallPotential
	m.FieldSwitchPosition = (1. - math.Pow(-m.Parameters.ConstEField*m.Parameters.CathodeFallLength/(m.Parameters.GetSheathFieldPower()*m.Parameters.CathodeFallPotential), 1./m.Parameters.GetSheathFieldPower())) * m.Parameters.CathodeFallLength
	m.cathodeRadius2 = m.Parameters.CathodeRadius * m.Parameters.CathodeRadius
	m.tubeRadius2 = m.Parameters.TubeRadius * m.Parameters.TubeRadius
	m.constFieldAcceleration = -constants.ElementaryCharge * m.Parameters.ConstEField / constants.ElectornMass

	m.sheathTimeOuterConstant = 2 * m.Parameters.CathodeFallLength * math.Sqrt(constants.ElectornMass/(2*constants.ElementaryCharge*m.Vs))
	m.sheathTimeInnerConstant = math.Sqrt(m.Parameters.CathodeFallPotential) / m.Parameters.CathodeFallLength
	m.constFieldTimeConstant = math.Sqrt(2*constants.ElectornMass/constants.ElementaryCharge) * (-m.Parameters.ConstEField)

	m.EFieldScale = -m.inverseCathodeFallLength * m.Parameters.CathodeFallPotential * (m.Parameters.GetSheathFieldPower() + 1.)
	m.FieldSwitchPotential = m.Parameters.ConstEField * (m.Parameters.SimulationLength() - m.FieldSwitchPosition)
	m.SheathPotentialCorrection = m.FieldSwitchPotential + m.Parameters.CathodeFallPotential*math.Pow((-m.Parameters.ConstEField*m.Parameters.CathodeFallLength/(m.Parameters.CathodeFallPotential*(m.Parameters.GetSheathFieldPower()+1))), (m.Parameters.GetSheathFieldPower()+1)/m.Parameters.GetSheathFieldPower())
}

func NewModel(parameters *config.ModelParameters) *Model {
	m := Model{}
	m.Parameters = *parameters
	m.InitVars()

	meanFreePath := 1. / (m.Parameters.CrossSectionsData().SurplusCrossSection() * m.Parameters.GasDensity)

	m.NumCells = int(5*m.Parameters.SimulationLength()/meanFreePath) + 1
	m.NumCellsMu = int(2./parameters.MuDiscretizationStep) + 1
	m.NumCellsE = int((-m.VfromL(0)+5)/m.Parameters.EnergyStep) + 1

	m.XStep = m.Parameters.SimulationLength() / float64(m.NumCells)

	m.TrajAveragedIonization = make([]utils.KahanSummable, m.NumCells+1)

	m.IonizationCSmulGasDensityLookupECell = make([]float64, m.NumCellsE+1)
	for p := range parameters.CrossSectionsData().Processes {
		if parameters.CrossSectionsData().Processes[p].Type == lxgata.IONIZATION {
			ionizationCS := &parameters.CrossSectionsData().Processes[p]
			for e := range m.IonizationCSmulGasDensityLookupECell {
				energy := (float64(e) + 0.5) * parameters.EnergyStep
				m.IonizationCSmulGasDensityLookupECell[e] = ionizationCS.CrossSectionAt(energy) * parameters.GasDensity
			}
			break
		}
	}

	m.CollisionAtCell = make(map[lxgata.CollisionType][]utils.Aggregation)
	m.IonizationsSumUpToCell = make([]utils.Aggregation, m.NumCells+1)
	for _, process := range parameters.CrossSectionsData().GetTypes() {
		m.CollisionAtCell[process] = make([]utils.Aggregation, m.NumCells+1)
	}
	m.DetailedCollisionAtCell = make(map[string][]utils.Aggregation) //utils.NewMutableMap[string, utils.AggregatedDistribution]() //make(map[string]utils.AggregatedDistribution)
	for p := range parameters.CrossSectionsData().Processes {
		for substring := range parameters.RequireCollisionRelativeMargin {
			collisionName := string(parameters.CrossSectionsData().Processes[p].Type) + parameters.CrossSectionsData().Processes[p].Outcome
			if strings.Contains(collisionName, substring) {
				m.DetailedCollisionAtCell[collisionName] = make([]utils.Aggregation, m.NumCells+1)
			}
		}
	}

	m.lookUpPotentialAtCellNode = make([]float64, m.NumCells+2)
	for cellNode := range m.lookUpPotentialAtCellNode {
		m.lookUpPotentialAtCellNode[cellNode] = m.VfromL(float64(cellNode) * m.XStep)
	}

	m.lookupCellNodeFromPotential = make([]int, m.NumCellsE+2)
	for e := range m.lookupCellNodeFromPotential {
		potentialEnergy := (float64(e) + 0.5) * parameters.EnergyStep
		m.lookupCellNodeFromPotential[e] = int(m.LfromV(-potentialEnergy) / m.XStep)
	}

	m.MeanFreePath = make([]utils.AggregationF, m.NumCells)

	numEnergyCells := m.NumCellsE + 10
	m.LookupEnergy = make([]float64, numEnergyCells)
	m.lookupMNumerator = make([]float64, numEnergyCells)
	for gridNode := range m.LookupEnergy {
		m.LookupEnergy[gridNode] = float64(gridNode) * m.Parameters.EnergyStep
		m.lookupMNumerator[gridNode] = m.Parameters.GasDensity * m.Parameters.CrossSectionsData().TotalCrossSectionAt(m.LookupEnergy[gridNode]) * math.Sqrt(m.LookupEnergy[gridNode])
	}

	m.MeanEnergy = make([]utils.AggregationF, m.NumCells+1)

	if m.Parameters.CalculateDistribution {
		m.DistributionXEMu = make([][][]float64, m.NumCells)
		for x := range m.DistributionXEMu {
			m.DistributionXEMu[x] = make([][]float64, m.NumCellsE)
			for e := range m.DistributionXEMu[x] {
				m.DistributionXEMu[x][e] = make([]float64, m.NumCellsMu)
			}
		}
	}

	m.DataHub = NewDataHub()

	return &m
}

func (m *Model) ReducedFieldAtCathode() float64 {
	return config.ReducedFieldAtCathode(m.Parameters.CathodeFallPotential, m.Parameters.CathodeFallLength, m.Parameters.GasDensity)
}

func (m *Model) ReducedFieldMidSheath() float64 {
	return config.ReducedFieldMidSheath(m.Parameters.CathodeFallPotential, m.Parameters.CathodeFallLength, m.Parameters.GasDensity)
}

type TrajectoryUpdate struct {
	startEnergy, endEnergy float64
	trajectory             TrajectoryConstants
	origin                 int
	weight                 int
}

type BoundaryConditionType int

const (
	None BoundaryConditionType = iota
	ArrivalAtCathode
	ArrivalAtSimEnd
	ArrivalAtTubeWall
)

type MisalignmentCause int

const (
	MCNone MisalignmentCause = iota
	SimEnd
	Cathode
	Turnover
)

func (m *Model) randomDistanceToWall() float64 {
	return 0
}

type PathRate struct {
	cell int
	rate float64
}

func (m *Model) nextCollision(p *Particle, tupd chan<- TrajectoryUpdate) (collisionDescription *lxgata.Collision, reversal bool, throwOut BoundaryConditionType, freePath utils.KahanSummable) {
	// timeToWall := m.randomDistanceToWall() / utils.EV2electronVelocity(p.trajectory.radialEnergy)
	R := utils.NewKahanSummable(rand.ExpFloat64()) //-math.Log(1. - rand.Float64()))

	if -p.getPotentialEnergy() < -m.Vs+m.SheathPotentialCorrection {
		throwOut = ArrivalAtCathode
		return
	}
	if -p.getPotentialEnergy() > 0 { //p.x > m.Parameters.SimulationLength() {
		throwOut = ArrivalAtSimEnd
		return
	}

	var currentCellIndex int
	var cachedVelocity float64
	var isVelocityCached = false

	alignedToEnergyGrid := false

	var pathStartEnergy, segmentEndEnergy float64
	if tupd != nil {
		pathStartEnergy = p.eKinetic
	}

	for {
		var lowEnergy, segmentStartEnergy, highEnergy float64
		var nextCellIndex, highEnergyCellIndex int
		var reversalOccured, arrivalAtCathode, arrivalAtSimEnd, arrivalAtTubeWall, highEnergyAligned bool
		fieldSwitch := false
		if p.mu < 0 {
			if !alignedToEnergyGrid {
				nextCellIndex = int(math.Floor(p.eKinetic / m.Parameters.EnergyStep))
				if nextCellIndex < m.NumCellsE {
					lowEnergy, highEnergy = m.LookupEnergy[nextCellIndex], p.eKinetic
				} else {
					lowEnergy, highEnergy = float64(nextCellIndex)*m.Parameters.EnergyStep, p.eKinetic
				}

			} else {
				if currentCellIndex == 0 {
					lowEnergy, highEnergy = -m.Parameters.EnergyStep, 0
				} else {
					if currentCellIndex < m.NumCellsE {
						lowEnergy, highEnergy = m.LookupEnergy[currentCellIndex-1], m.LookupEnergy[currentCellIndex]
					} else {
						lowEnergy, highEnergy = float64(currentCellIndex-1)*m.Parameters.EnergyStep, float64(currentCellIndex)*m.Parameters.EnergyStep
					}
				}

				highEnergyCellIndex, highEnergyAligned = currentCellIndex, true
				nextCellIndex = currentCellIndex - 1
			}

			//check for reversal
			if lowEnergy < p.trajectory.radialEnergy {
				lowEnergy, alignedToEnergyGrid = p.trajectory.radialEnergy, false
				nextCellIndex = currentCellIndex
				reversalOccured = true
			}

			//check for cathode return
			if lowEnergy-p.trajectory.totEnergy < m.SheathPotentialCorrection-m.Vs {
				arrivalAtCathode = true
				lowEnergy = p.trajectory.totEnergy - m.Vs + m.SheathPotentialCorrection
			}

			// fieldSwitchEnergy := p.trajectory.totEnergy + m.FieldSwitchPotential
			// if lowEnergy < fieldSwitchEnergy && fieldSwitchEnergy < highEnergy {
			// 	lowEnergy = fieldSwitchEnergy - 1e-9
			// 	nextCellIndex = currentCellIndex
			// 	alignedToEnergyGrid = false
			// 	fieldSwitch = true
			// }

			segmentStartEnergy, segmentEndEnergy = highEnergy, lowEnergy
		} else {
			if !alignedToEnergyGrid {
				nextCellIndex = int(p.eKinetic/m.Parameters.EnergyStep) + 1
				if nextCellIndex < m.NumCellsE {
					lowEnergy, highEnergy = p.eKinetic, m.LookupEnergy[nextCellIndex]
				} else {
					lowEnergy, highEnergy = p.eKinetic, float64(nextCellIndex)*m.Parameters.EnergyStep
				}

			} else {
				if currentCellIndex < m.NumCellsE {
					lowEnergy, highEnergy = m.LookupEnergy[currentCellIndex], m.LookupEnergy[currentCellIndex+1]
				} else {
					lowEnergy, highEnergy = float64(currentCellIndex)*m.Parameters.EnergyStep, float64(currentCellIndex+1)*m.Parameters.EnergyStep
				}

				nextCellIndex = currentCellIndex + 1
			}

			//check for gap end arrival
			if p.trajectory.totEnergy <= highEnergy {
				arrivalAtSimEnd = true
				highEnergy, highEnergyAligned = p.trajectory.totEnergy, false
			} else {
				highEnergyCellIndex, highEnergyAligned = nextCellIndex, true
			}
			// fieldSwitchEnergy := p.trajectory.totEnergy + m.FieldSwitchPotential
			// if lowEnergy < fieldSwitchEnergy && fieldSwitchEnergy < highEnergy {
			// 	highEnergy = fieldSwitchEnergy + 1e-9
			// 	highEnergyAligned = false
			// 	nextCellIndex = currentCellIndex
			// 	alignedToEnergyGrid = false
			// 	fieldSwitch = true
			// }

			segmentStartEnergy, segmentEndEnergy = lowEnergy, highEnergy
		}

		var M, G float64
		if highEnergyAligned && highEnergyCellIndex < m.NumCellsE {
			M = m.lookupMNumerator[highEnergyCellIndex] / -m.EFieldFromPotential(-(p.trajectory.totEnergy - highEnergy))
		} else {
			M = m.Parameters.GasDensity * m.Parameters.CrossSectionsData().TotalCrossSectionAt(highEnergy) * math.Sqrt(highEnergy) / -m.EFieldFromPotential(-(p.trajectory.totEnergy - highEnergy))
		}

		var higherVelocity, lowerVelocity float64
		if isVelocityCached && !fieldSwitch {
			if p.mu < 0 {
				higherVelocity, lowerVelocity = cachedVelocity, math.Sqrt(lowEnergy-p.trajectory.radialEnergy)
			} else {
				higherVelocity, lowerVelocity = math.Sqrt(highEnergy-p.trajectory.radialEnergy), cachedVelocity
			}
		} else {
			higherVelocity, lowerVelocity = math.Sqrt(highEnergy-p.trajectory.radialEnergy), math.Sqrt(lowEnergy-p.trajectory.radialEnergy)
		}
		if p.mu < 0 {
			cachedVelocity = lowerVelocity
		} else {
			cachedVelocity = higherVelocity
		}
		isVelocityCached = true

		G = 2 * M * (higherVelocity - lowerVelocity)

		collisionOccured := false
		if G < R.Val() {
			R.Add(-G)
			if m.Parameters.MeanFreePath {
				meanEnergy := 0.5 * (segmentStartEnergy + segmentEndEnergy)
				meanMu := math.Sqrt((meanEnergy - p.trajectory.radialEnergy) / meanEnergy)
				deltaX := m.LfromV(highEnergy-p.trajectory.totEnergy) - m.LfromV(lowEnergy-p.trajectory.totEnergy)
				deltaPath := deltaX / meanMu
				if !math.IsNaN(deltaPath) && !math.IsInf(deltaPath, 0) {
					freePath.Add(deltaX / meanMu)
				}

			}
		} else {
			{
				var delta float64
				if p.mu < 0 {
					delta = math.Sqrt(segmentStartEnergy-p.trajectory.radialEnergy) - R.Val()/(2.*M) // == sqrt(pColl - p.trajectory.eStar) i.e. abs(speed) equivalent along x axis
					if delta < 0 {
						R.Add(-2 * M * math.Sqrt(segmentStartEnergy-p.trajectory.radialEnergy))
						segmentEndEnergy = p.trajectory.radialEnergy
						reversalOccured = true
					}
					if math.IsNaN(delta) {
						panic("delta is NaN")
					}
				} else {
					delta = R.Val()/(2.*M) + math.Sqrt(segmentStartEnergy-p.trajectory.radialEnergy)
					if math.IsNaN(delta) {
						panic("delta is NaN")
					}
				}
				if !reversalOccured {
					collisionEnergy := math.FMA(delta, delta, p.trajectory.radialEnergy) // R = 2M[sqrt(p.e-p.trajectory.eStar) - sqrt(eColl - p.trajectory.eStar)]
					segmentEndEnergy = collisionEnergy
					p.setEnergy(collisionEnergy, m, true)
					if p.trajectory.totEnergy < collisionEnergy || p.getPotentialEnergy() < 0 {
						arrivalAtSimEnd = true
						segmentEndEnergy = p.trajectory.totEnergy
					} else if -p.getPotentialEnergy() < -m.Vs+m.SheathPotentialCorrection {
						arrivalAtCathode = true
						segmentEndEnergy = p.trajectory.totEnergy - m.Vs + m.SheathPotentialCorrection
					} else {
						if !m.Parameters.Volumetric { //|| p.y*p.y+p.z*p.z < m.tubeRadius2 {
							var totalCrossSectionPrimed = M * math.Abs(m.EFieldFromPotential(collisionEnergy-p.trajectory.totEnergy)) / (math.Sqrt(collisionEnergy) * m.Parameters.GasDensity)
							// collisionDescription = m.collisionSelector(collisionEnergy, p.x, M)
							collisionDescription = m.Parameters.CrossSectionsData().SampleWithNullCollision(collisionEnergy, totalCrossSectionPrimed)
							collisionOccured = true
						} else {
							arrivalAtTubeWall = true
						}
					}
				}
				if m.Parameters.MeanFreePath {
					if p.mu < 0 {
						meanEnergy := 0.5 * (segmentStartEnergy + segmentEndEnergy)
						meanMu := math.Sqrt((meanEnergy - p.trajectory.radialEnergy) / meanEnergy)
						deltaX := m.LfromV(highEnergy-p.trajectory.totEnergy) - m.LfromV(segmentEndEnergy-p.trajectory.totEnergy)
						deltaPath := deltaX / meanMu
						if !math.IsNaN(deltaPath) && !math.IsInf(deltaPath, 0) {
							freePath.Add(deltaX / meanMu)
						}
					} else {
						meanEnergy := 0.5 * (segmentStartEnergy + segmentEndEnergy)
						meanMu := math.Sqrt((meanEnergy - p.trajectory.radialEnergy) / meanEnergy)
						deltaX := m.LfromV(segmentEndEnergy-p.trajectory.totEnergy) - m.LfromV(lowEnergy-p.trajectory.totEnergy)
						deltaPath := deltaX / meanMu
						if !math.IsNaN(deltaPath) && !math.IsInf(deltaPath, 0) {
							freePath.Add(deltaX / meanMu)
						}
					}
				}
			}
		}

		if collisionOccured || arrivalAtTubeWall {
			break
		}

		if arrivalAtCathode {
			throwOut = ArrivalAtCathode
			break
		}

		if arrivalAtSimEnd {
			if m.Parameters.ParallelPlaneHollowCathode {
				if p.mu < 0 {
					println("DEBUG: arrival at half-gap from beyond")
				}
				p.mu = -p.mu
				// p.eKinetic = p.trajectory.totEnergy
				p.setEnergy(p.trajectory.totEnergy, m, true)
				alignedToEnergyGrid = false
			} else {
				if m.Parameters.AnodeBackscatteringCoefficient0 != 0 {
					backscatteringCoefficient := m.Parameters.AnodeBackscatteringCoefficient0 * math.Exp(m.Parameters.AnodeBackscatteringCoefficientB*(1-p.mu))

					p.setEnergy(p.trajectory.totEnergy, m, true)
					if rand.Float64() < backscatteringCoefficient { //&& p.y*p.y+p.z*p.z < m.cathodeRadius2 {
						p.mu = -p.mu
						alignedToEnergyGrid = false
						if m.Parameters.AnodeBackscatteringEnergyLossFraction != 0 {
							p.eKinetic *= (1. - m.Parameters.AnodeBackscatteringEnergyLossFraction)
						} else if p.eKinetic > 9 && rand.Float64() > 0.5 {
							p.eKinetic = p.eKinetic - 9
						}
						p.recalcParams(m)
					} else {
						throwOut = ArrivalAtSimEnd
						break
					}
				} else {
					throwOut = ArrivalAtSimEnd
					p.eKinetic = p.trajectory.totEnergy
					p.mu = math.Sqrt(p.getAxialEnergy() / p.trajectory.totEnergy)
					break
				}

			}
		} else if reversalOccured {
			isVelocityCached = false
			p.mu = +0.
			p.eKinetic = p.trajectory.radialEnergy
			alignedToEnergyGrid = false
			// pathStartEnergy = p.trajectory.radialEnergy
			reversal = true
		} else if !fieldSwitch {
			alignedToEnergyGrid = true
		}
		if m.Parameters.Volumetric && !math.IsInf(m.tubeRadius2, 0) { //&& currentCellIndex == wallCollisionEnergyStep
			throwOut = ArrivalAtTubeWall
			break
		}
		if !alignedToEnergyGrid {
			p.setEnergy(segmentEndEnergy, m, true)
		}

		currentCellIndex = nextCellIndex
	}
	if tupd != nil {
		if reversal {
			tupd <- TrajectoryUpdate{
				startEnergy: pathStartEnergy,
				endEnergy:   p.trajectory.radialEnergy,
				trajectory:  p.trajectory,
				origin:      p.origin,
				weight:      p.weight,
			}
			tupd <- TrajectoryUpdate{
				startEnergy: p.trajectory.radialEnergy,
				endEnergy:   segmentEndEnergy,
				trajectory:  p.trajectory,
				origin:      p.origin,
				weight:      p.weight,
			}
		} else {
			tupd <- TrajectoryUpdate{
				startEnergy: pathStartEnergy,
				endEnergy:   segmentEndEnergy,
				trajectory:  p.trajectory,
				origin:      p.origin,
				weight:      p.weight,
			}
		}
	}
	return
}

type CollisionEvent struct {
	x          int
	energyLoss float64
	collType   lxgata.CollisionType
	outcome    string
	origin     int
	weight     int
}

type MeanFreePathUpdate struct {
	freePath                     utils.KahanSummable
	fillStartIndex, fillEndIndex int
}

type AnodeArrival struct {
	weight int
}

type CathodeArrival struct {
	weight int
	origin int
}

type DriftSegment struct {
	isAttachment bool
	origin       int
	dx           float64
}

type IonizingElectron struct {
	origin int
	weight int
}

type CountEvent struct { //ionization or attachment
	isAttachment bool
	origin       int
	weight       int
}

func (m *Model) Run(particlesToLaunch func(*Model) int, electrons bool, logger messages.Logger) {
	CollisionAtCellPerAvalanche := make(map[lxgata.CollisionType][][]int)
	for process := range m.CollisionAtCell {
		CollisionAtCellPerAvalanche[process] = make([][]int, m.NumCells+1)
	}
	DetailedCollisionAtCell := make(map[string][][]int)
	for process := range m.DetailedCollisionAtCell {
		DetailedCollisionAtCell[process] = make([][]int, m.NumCells+1)
	}

	var energyDeposition float64
	var numberOfEDEvents int

	nParticles := particlesToLaunch(m)
	for nParticles > 0 {
		m.PurgeMetrics()
		cfCap := nParticles * 1000
		computeFlow := make(chan *Particle, cfCap)
		var computeWg, stateWg sync.WaitGroup

		if electrons {

		}

		ionizationCounters := make([]int, nParticles)
		attachmentCounters := make([]int, nParticles)

		collFlow := make(chan CollisionEvent, 1000*m.Parameters.NParticles)
		stateWg.Add(1)
		go func() {
			for collision := range collFlow {
				if 0 < collision.x && collision.x < m.NumCells {
					CollisionAtCellPerAvalanche[collision.collType][collision.x][collision.origin-m.TotalElectronsEmittedOnCathode] += collision.weight
					collName := string(collision.collType) + collision.outcome
					if _, exist := DetailedCollisionAtCell[collName]; exist {
						DetailedCollisionAtCell[collName][collision.x][collision.origin-m.TotalElectronsEmittedOnCathode] += collision.weight
					}
				}
			}
			stateWg.Done()
		}()

		var trajFlow chan TrajectoryUpdate
		if m.Parameters.TrajectoryAveragedRates || m.Parameters.CalculateDistribution {
			// for x := range TimeAtCellPerAvalanche {
			// 	TimeAtCellPerAvalanche[x] = make([]utils.KahanSummable, nElectrons)
			// }

			trajFlow = make(chan TrajectoryUpdate, 100000)
			stateWg.Add(1)
			go func() {
				// if m.Parameters.TrajectoryAveragedIonization {
				// 	for tUpdate := range trajFlow {
				// 		es, ee := min(tUpdate.startEnergy, tUpdate.endEnergy), max(tUpdate.startEnergy, tUpdate.endEnergy)
				// 		m.RateIncrement(es, float64(int(math.Ceil(es/m.Parameters.EnergyStep)))*m.Parameters.EnergyStep, tUpdate.trajectory.totEnergy, m.TrajAveragedIonization)
				// 		topEi := int(math.Floor(ee / m.Parameters.EnergyStep))
				// 		for i := int(math.Ceil(es / m.Parameters.EnergyStep)); i < topEi; i++ {
				// 			m.RateIncrement(float64(i)*m.Parameters.EnergyStep, float64(i+1)*m.Parameters.EnergyStep, tUpdate.trajectory.totEnergy, m.TrajAveragedIonization)
				// 		}
				// 	}
				// }
				if m.Parameters.TrajectoryAveragedRates || m.Parameters.CalculateDistribution {
					for tUpdate := range trajFlow {
						startEnergy, endEnergy := min(tUpdate.startEnergy, tUpdate.endEnergy), max(tUpdate.startEnergy, tUpdate.endEnergy)
						if m.Parameters.TrajectoryAveragedRates && endEnergy > m.Parameters.LowerEnergyThreshold() {
							eCell := int(math.Ceil(max(m.Parameters.LowerEnergyThreshold(), startEnergy) / m.Parameters.EnergyStep))
							endCellIndex := int(math.Floor(endEnergy/m.Parameters.EnergyStep)) - 1
							for ; eCell < endCellIndex; eCell++ {
								eMid := m.Parameters.EnergyStep * (float64(eCell) + 0.5)
								invMu := math.Sqrt(eMid / (eMid - tUpdate.trajectory.radialEnergy))
								potentialMid := eMid - tUpdate.trajectory.totEnergy
								field := m.EFieldFromPotential(potentialMid)
								ncs := m.IonizationCSmulGasDensityLookupECell[eCell]
								for range tUpdate.weight {
									m.TrajAveragedIonization[int(m.LfromV(potentialMid)/m.XStep)].Add(ncs * m.Parameters.EnergyStep * invMu / -field)
								}
							}
							{
								//interval boundaries
								ingressCellEnergy := math.Ceil(startEnergy/m.Parameters.EnergyStep) * m.Parameters.EnergyStep
								egressCellEnergy := math.Floor(endEnergy/m.Parameters.EnergyStep) * m.Parameters.EnergyStep
								if ingressCellEnergy <= endEnergy {
									startStep := ingressCellEnergy - startEnergy
									eMid := 0.5 * (startEnergy + ingressCellEnergy)
									invMu := math.Sqrt(eMid / (eMid - tUpdate.trajectory.radialEnergy))
									potentialMid := eMid - tUpdate.trajectory.totEnergy
									x := m.LfromV(potentialMid)
									if x > m.Parameters.CathodeFallLength {
										invMu = math.Pi * 0.5
									}
									field := m.EFieldFromPotential(potentialMid)
									ncs := m.IonizationCSmulGasDensityLookupECell[eCell]
									for range tUpdate.weight {
										m.TrajAveragedIonization[int(m.LfromV(potentialMid)/m.XStep)].Add(ncs * startStep * invMu / -field)
									}
								}
								if egressCellEnergy >= startEnergy {
									endStep := endEnergy - egressCellEnergy
									eMid := 0.5 * (egressCellEnergy + endEnergy)
									invMu := math.Sqrt(eMid / (eMid - tUpdate.trajectory.radialEnergy))
									potentialMid := eMid - tUpdate.trajectory.totEnergy
									x := m.LfromV(potentialMid)
									if x > m.Parameters.CathodeFallLength {
										invMu = math.Pi * 0.5
									}
									field := m.EFieldFromPotential(potentialMid)
									ncs := m.IonizationCSmulGasDensityLookupECell[eCell]
									for range tUpdate.weight {
										m.TrajAveragedIonization[int(m.LfromV(potentialMid)/m.XStep)].Add(ncs * endStep * invMu / -field)
									}
								}
								if ingressCellEnergy > endEnergy && egressCellEnergy < startEnergy {
									interCellStep := endEnergy - startEnergy
									eMid := 0.5 * (startEnergy + endEnergy)
									invMu := math.Sqrt(eMid / (eMid - tUpdate.trajectory.radialEnergy))
									potentialMid := eMid - tUpdate.trajectory.totEnergy
									x := m.LfromV(potentialMid)
									if x > m.Parameters.CathodeFallLength {
										invMu = math.Pi * 0.5
									}
									field := m.EFieldFromPotential(potentialMid)
									ncs := m.IonizationCSmulGasDensityLookupECell[eCell]
									for range tUpdate.weight {
										m.TrajAveragedIonization[int(x/m.XStep)].Add(ncs * interCellStep * invMu / -field)
									}
								}
							}
						}
						if m.Parameters.CalculateDistribution {
							xStart := m.LfromV(tUpdate.startEnergy - tUpdate.trajectory.totEnergy)
							xEnd := m.LfromV(tUpdate.endEnergy - tUpdate.trajectory.totEnergy)
							xStart, xEnd = min(xStart, xEnd), max(xStart, xEnd)
							xStartIndex := int(math.Ceil(xStart / m.XStep))
							xEndIndex := int(math.Floor(xEnd / m.XStep))
							for xIndex := xStartIndex; xIndex <= xEndIndex; xIndex++ {
								energy := tUpdate.trajectory.totEnergy + m.lookUpPotentialAtCellNode[xIndex]
								eIndex := int(energy / m.Parameters.EnergyStep)
								mu := math.Copysign(math.Sqrt((energy-tUpdate.trajectory.radialEnergy)/energy), tUpdate.endEnergy-tUpdate.startEnergy)
								muIndex := int((mu + 1 + 0.5*m.Parameters.MuDiscretizationStep) / m.Parameters.MuDiscretizationStep)
								if xIndex < m.NumCells && eIndex < m.NumCellsE && muIndex < m.NumCellsMu {
									for range tUpdate.weight {
										m.MeanEnergy[xIndex].AppendF(energy)
									}
									m.DistributionXEMu[xIndex][eIndex][muIndex] += float64(tUpdate.weight)
								}
							}
						}
					}
				}
				stateWg.Done()
			}()
		}

		for _, process := range m.Parameters.CrossSectionsData().GetTypes() {
			for x := range CollisionAtCellPerAvalanche[process] {
				CollisionAtCellPerAvalanche[process][x] = make([]int, nParticles)
			}
		}

		var freePathFlow chan MeanFreePathUpdate
		if m.Parameters.MeanFreePath {
			freePathFlow = make(chan MeanFreePathUpdate)
			stateWg.Add(1)
			go func() {
				for fpu := range freePathFlow {
					for i := fpu.fillStartIndex; i <= fpu.fillEndIndex && i < m.NumCells; i++ {
						m.MeanFreePath[i].Append(fpu.freePath)
					}
					m.GlobalMeanFreePath.Append(fpu.freePath)
				}
				stateWg.Done()
			}()
		}

		electronsReturned := make([]int, nParticles)
		erFlow := make(chan CathodeArrival)
		stateWg.Add(1)
		go func() {
			for cathodeReturn := range erFlow {
				electronsReturned[cathodeReturn.origin-m.TotalElectronsEmittedOnCathode] += cathodeReturn.weight
			}
			stateWg.Done()
		}()
		ionizingElectrons := make([]int, nParticles)
		izFlow := make(chan IonizingElectron)
		stateWg.Add(1)
		go func() {
			for ie := range izFlow {
				ionizingElectrons[ie.origin-m.TotalElectronsEmittedOnCathode] += ie.weight
			}
			stateWg.Done()
		}()

		anodeFlow := make(chan AnodeArrival)
		stateWg.Add(1)
		go func() {
			for anodeArrival := range anodeFlow {
				m.AnodeElectronCounter += anodeArrival.weight
			}
			stateWg.Done()
		}()

		countFlow := make(chan CountEvent)
		stateWg.Add(1)
		go func() {
			for event := range countFlow {
				if event.isAttachment {
					attachmentCounters[event.origin-m.TotalElectronsEmittedOnCathode] += event.weight
				} else {
					ionizationCounters[event.origin-m.TotalElectronsEmittedOnCathode] += event.weight
				}
			}
			stateWg.Done()
		}()

		for origin := range nParticles {
			particle := m.newParticle(origin + m.TotalElectronsEmittedOnCathode)
			if m.Parameters.CalculateDistribution {
				m.DistributionXEMu[0][int(particle.eKinetic/m.Parameters.EnergyStep)][int((particle.mu+1)/m.Parameters.MuDiscretizationStep)]++
			}
			computeWg.Add(1)
			computeFlow <- &particle
		}

		// status := []string{"//", "==", "\\\\", "||"}
		closeCallback := func() {}
		if !m.Parameters.SupressSpinner {
			closeCallback = logger.Busy("Avalanches")
		}

		reweightRequest := false

		go func() {
			computeWg.Wait()

			close(computeFlow)
			close(collFlow)
			if m.Parameters.TrajectoryAveragedRates || m.Parameters.CalculateDistribution {
				close(trajFlow)
			}

			if m.Parameters.MeanFreePath {
				close(freePathFlow)
			}
			close(erFlow)
			close(izFlow)
			close(anodeFlow)
			close(countFlow)
		}()

		for {
			var workerWG sync.WaitGroup
			for range m.Parameters.Threads() - 1 {
				workerWG.Add(1)
				go func() {
					for particlePtr := range computeFlow {
						lowerEnergyThreshold := m.Parameters.LowerEnergyThreshold()
						attachmentLoss := false
						for (lowerEnergyThreshold < particlePtr.trajectory.totEnergy || m.Parameters.CalculateDistribution) &&
							(!m.Parameters.ReturningElectrons || particlePtr.trajectory.totEnergy+m.VfromL(0) > 0) {
							if len(computeFlow)*10/(cap(computeFlow)*8) >= 1 || reweightRequest {
								computeFlow <- particlePtr
								computeWg.Add(1)
								reweightRequest = true
								break
							}
							var fillStartIndex, fillEndIndex int
							if m.Parameters.MeanFreePath {
								fillStartIndex = m.LNodefromVCached(-particlePtr.getPotentialEnergy())
							}
							collision, reversal, throwOut, freePath := m.nextCollision(particlePtr, trajFlow /*, stateflow*/)
							if m.Parameters.MeanFreePath {
								fillEndIndex = m.LNodefromVCached(-particlePtr.getPotentialEnergy())
								fillStartIndex, fillEndIndex = min(fillStartIndex, fillEndIndex), max(fillStartIndex, fillEndIndex)
								if reversal {
									fillStartIndex = int(particlePtr.trajectory.getTurnaroundX(m) / m.XStep)
								}
								freePathFlow <- MeanFreePathUpdate{
									freePath:       freePath,
									fillStartIndex: fillStartIndex,
									fillEndIndex:   fillEndIndex,
								}
							}

							if throwOut != None {
								switch throwOut {
								case ArrivalAtCathode:
									erFlow <- CathodeArrival{weight: particlePtr.weight, origin: particlePtr.origin}
								case ArrivalAtSimEnd:
									anodeFlow <- AnodeArrival{weight: particlePtr.weight}
								}
								break
							}
							if collision != nil {
								energyLoss := collision.Threshold
								cosChiScattered := m.Parameters.CrossSectionsData().SampleScatteringAngleCos(particlePtr.eKinetic, energyLoss, collision.Type, lxgata.IgnoreAtomicNumber, collision.Species)
								phi := 2. * math.Pi * rand.Float64()

								particlePtr.eKinetic -= energyLoss
								switch collision.Type {
								case lxgata.ELASTIC:
									energyLoss = particlePtr.eKinetic * 2. * collision.MassRatio * (1. - cosChiScattered)
									particlePtr.eKinetic -= energyLoss

								case lxgata.IONIZATION:
									particlePtr.generation += 1
									ejected := *particlePtr
									ejected.ejectedFromIonization = true
									particlePtr.producedIonization = true

									availableEnergy := particlePtr.eKinetic

									switch m.Parameters.GetIonizationEnergySharingMode() {
									case config.Equal:
										ejected.eKinetic = 0.5 * particlePtr.eKinetic
									case config.Opal:
										w := utils.OpalIonizationShapeParameters[collision.Species]
										ejected.eKinetic = w * math.Tan(rand.Float64()*math.Atan(particlePtr.eKinetic/(2.*w)))
									case config.UniformRandom:
										ejected.eKinetic = particlePtr.eKinetic * rand.Float64()
									case config.Zero:
										ejected.eKinetic = 0
									default:
										panic("unexpected config.IonizationEnergySharing")
									}
									particlePtr.eKinetic = availableEnergy - ejected.eKinetic

									var cosChiEjected float64

									switch m.Parameters.GetIonizationScatteringMode() {
									case config.Boeuf:
										cosChiScattered = math.Sqrt(particlePtr.eKinetic / availableEnergy)
										cosChiEjected = math.Sqrt(ejected.eKinetic / availableEnergy)
									case config.Born:
										cosChiScattered = lxgata.BornScatteringAngleSample(particlePtr.eKinetic+energyLoss, energyLoss)
										cosChiEjected = lxgata.BornScatteringAngleSample(ejected.eKinetic+energyLoss, energyLoss)
									case config.BornFullLoss:
										cosChiScattered = lxgata.BornScatteringAngleSample(availableEnergy+energyLoss, energyLoss+ejected.eKinetic)
										cosChiEjected = lxgata.BornScatteringAngleSample(availableEnergy+energyLoss, energyLoss+particlePtr.eKinetic)
									case config.Coulomb:
										cosChiScattered = lxgata.CoulombScatteringAngleSample(particlePtr.eKinetic+energyLoss, m.Parameters.CrossSectionsData().Species[collision.Species].UParameter, energyLoss, lxgata.IgnoreAtomicNumber)
										cosChiEjected = 1. - 2.*rand.Float64()
									case config.CoulombM:
										cosChiScattered = lxgata.CoulombScatteringAngleSample(particlePtr.eKinetic+energyLoss, m.Parameters.CrossSectionsData().Species[collision.Species].UParameter, energyLoss, lxgata.IgnoreAtomicNumber)
										cosChiEjected = lxgata.CoulombScatteringAngleSample(ejected.eKinetic+energyLoss, m.Parameters.CrossSectionsData().Species[collision.Species].UParameter, energyLoss, lxgata.IgnoreAtomicNumber)
									case config.Isotropic:
										cosChiScattered = 1. - 2.*rand.Float64()
										cosChiEjected = 1. - 2.*rand.Float64()
									case config.Mixed:
										eScattered, eEjected := max(particlePtr.eKinetic, ejected.eKinetic), min(particlePtr.eKinetic, ejected.eKinetic)
										particlePtr.eKinetic, ejected.eKinetic = eScattered, eEjected
										cosChiScattered = lxgata.BornScatteringAngleSample(particlePtr.eKinetic+energyLoss, energyLoss)
										cosChiEjected = 1. - 2.*rand.Float64()
									case config.MixedFullLoss:
										eScattered, eEjected := max(particlePtr.eKinetic, ejected.eKinetic), min(particlePtr.eKinetic, ejected.eKinetic)
										particlePtr.eKinetic, ejected.eKinetic = eScattered, eEjected
										cosChiScattered = lxgata.BornScatteringAngleSample(availableEnergy+energyLoss, energyLoss+ejected.eKinetic)
										cosChiEjected = 1. - 2.*rand.Float64()
									default:
										panic("unexpected config.IonizationScattering")
									}

									countFlow <- CountEvent{isAttachment: false, origin: particlePtr.origin, weight: particlePtr.weight}

									ejected.redirect(cosChiEjected, math.Cos(phi+math.Pi), m)
									if ejected.trajectory.totEnergy > lowerEnergyThreshold {
										computeWg.Add(1)
										computeFlow <- &ejected
									} else {
										anodeFlow <- AnodeArrival{weight: ejected.weight}
									}

								case lxgata.ATTACHMENT:
									countFlow <- CountEvent{isAttachment: true, origin: particlePtr.origin, weight: particlePtr.weight}
								case lxgata.EFFECTIVE:
								case lxgata.EXCITATION: //energy is lost with threshold decrement
								case lxgata.ROTATION:
								case lxgata.DEEXCITATION: //energy is gained with (negative) threshold decrement

								default:
									panic(fmt.Sprintf("unexpected lxgata.CollisionType: %#v", collision.Type))
								}
								// if !m.Parameters.Volumetric || (particlePtr.y*particlePtr.y+particlePtr.z*particlePtr.z < m.cathodeRadius2) || particlePtr.x < m.Parameters.CathodeFallLength {
								if _, count := m.Parameters.GetCollisionTypesToStore()[collision.Type]; count {
									collFlow <- CollisionEvent{
										x:          m.LNodefromVCached(particlePtr.potential),
										energyLoss: energyLoss,
										collType:   collision.Type,
										origin:     particlePtr.origin,
										outcome:    collision.Outcome,
										weight:     particlePtr.weight,
									}
								}

								// } // else {
								// 	outcollFlow <- CollisionEvent{
								// 		x:          int(particlePtr.x / m.XStep),
								// 		energyLoss: energyLoss,
								// 		collType:   collision.Type,
								// 		origin:     particlePtr.origin,
								// 		outcome:    collision.Outcome,
								// 	}
								// }

								if collision.Type == lxgata.ATTACHMENT {
									attachmentLoss = true
									break
								}
								particlePtr.redirect(cosChiScattered, math.Cos(phi), m)
							}
						}

						if particlePtr.trajectory.totEnergy < lowerEnergyThreshold && !attachmentLoss {
							anodeFlow <- AnodeArrival{weight: particlePtr.weight}
						}

						if !particlePtr.ejectedFromIonization && particlePtr.producedIonization {
							izFlow <- IonizingElectron{origin: particlePtr.origin, weight: particlePtr.weight}
						}
						computeWg.Done()
						if reweightRequest {
							break
						}
					}
					workerWG.Done()
				}()
			}
			workerWG.Wait()
			if reweightRequest {
				newCF := make(chan *Particle, cfCap)
				for particlePtr := range computeFlow {
					if rand.Float64() < 0.5 {
						particlePtr.weight *= 2
						newCF <- particlePtr
					} else {
						computeWg.Done()
					}
					if len(computeFlow) == 0 {
						break
					}
				}
				close(computeFlow)
				computeFlow = newCF
				reweightRequest = false
			} else {
				break
			}
		}

		stateWg.Wait()

		for process := range m.CollisionAtCell {
			for x := range m.NumCells {
				m.CollisionAtCell[process][x].Update(CollisionAtCellPerAvalanche[process][x])
			}
		}
		for process := range m.DetailedCollisionAtCell {
			for x := range m.NumCells {
				m.DetailedCollisionAtCell[process][x].Update(DetailedCollisionAtCell[process][x])
			}
		}

		PerAvalancheCollisionSumsUpToCell := make([]int, nParticles)
		for x := range m.IonizationsSumUpToCell {
			for electron := range nParticles {
				PerAvalancheCollisionSumsUpToCell[electron] += CollisionAtCellPerAvalanche[lxgata.IONIZATION][x][electron]
			}
			m.IonizationsSumUpToCell[x].Update(PerAvalancheCollisionSumsUpToCell)

		}
		m.CathodeElectronCounter += utils.SumIntSlice(electronsReturned)
		m.ElectronsReturned.Update(electronsReturned)
		m.IonizingCathodeElectrons.Update(ionizingElectrons)
		m.TotalElectronsEmittedOnCathode += nParticles
		m.IonizationCounters = append(m.IonizationCounters, ionizationCounters...)
		m.AttachmentCounters = append(m.AttachmentCounters, attachmentCounters...)
		nParticles = particlesToLaunch(m)
		closeCallback()
	}

	if numberOfEDEvents > 0 {
		m.MeanElectronEnergyAtAnode = energyDeposition / float64(numberOfEDEvents)
	}
	print("\r")
}
