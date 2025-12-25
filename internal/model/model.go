package model

import (
	"fmt"
	"math"
	"math/rand"
	"strings"
	"sync"

	"github.com/wildstyl3r/lxgata"
	"github.com/wildstyl3r/stmc/internal/config"
	"github.com/wildstyl3r/stmc/internal/constants"
	"github.com/wildstyl3r/stmc/internal/utils"
)

type Model struct {
	Parameters               config.ModelParameters
	Vc                       float64 // cathode fall potential
	EFieldA, EFieldBb        float64
	motionK, motionW         float64
	VOffset                  float64
	inverseConstEField       float64
	inverseCathodeFallLength float64
	sheathTimeOuterConstant  float64
	sheathTimeInnerConstant  float64
	constFieldTimeConstant   float64
	constFieldAcceleration   float64

	cathodeRadius2, tubeRadius2 float64

	Va float64 // additional voltage to avoid numeric negative energy beyond cathode fall

	NumCells, NumCellsMu, NumCellsE int

	XStep float64

	lookUpPotentialAtCellNode []float64
	LookupEnergy              []float64
	lookupMNumerator          []float64

	DistributionXEMu        [][][]float64
	RateLookups             map[string][]float64
	ionizationCS            *lxgata.Collision
	TimeAtCell              []utils.AggregationF
	MeanFreePath            []utils.AggregationF
	GlobalMeanFreePath      utils.AggregationF
	CollisionAtCell         map[lxgata.CollisionType][]utils.Aggregation
	DetailedCollisionAtCell map[string][]utils.Aggregation
	MeanEnergy              []utils.AggregationF
	IonizationsSumUpToCell  []utils.Aggregation

	TotalElectronsEmittedOnCathode int
	ElectronsReturned              utils.Aggregation

	MeanElectronEnergyAtAnode float64

	DataHub *DataHubType
}

func (m *Model) CoreResult() utils.CoreResult {
	electronsReturned, electronsReturnedMargin := m.ElectronsReturned.MeanWithErrorMargin(0.95)
	return utils.CoreResult{
		ModelName:                  m.Parameters.PrototypeName(),
		ReducedFieldAtCathode:      m.ReducedFieldAtCathode(),
		ReducedFieldAtSheathCenter: m.ReducedFieldMidSheath(),
		CathodeFallLength:          m.Parameters.CathodeFallLength,
		PressureCathodeFallLength:  m.Parameters.CathodeFallLength * m.Parameters.Pressure,
		GlobalMeanFreePath:         m.GlobalMeanFreePath.Mean(),
		Voltage:                    m.Parameters.CathodeFallPotential,
		Pressure:                   m.Parameters.Pressure,
		ElectronsReturned:          electronsReturned,
		ElectronsReturnedMargin:    electronsReturnedMargin,
	}
}

func (m *Model) InitVars() {
	m.Vc = math.Abs(m.Parameters.CathodeFallPotential)
	m.Va = math.Abs(m.Parameters.ConstEField * (m.Parameters.SimulationLength() - m.Parameters.CathodeFallLength))
	m.inverseConstEField = 1. / m.Parameters.ConstEField
	m.inverseCathodeFallLength = 1. / m.Parameters.CathodeFallLength
	m.cathodeRadius2 = m.Parameters.CathodeRadius * m.Parameters.CathodeRadius
	m.tubeRadius2 = m.Parameters.TubeRadius * m.Parameters.TubeRadius
	m.constFieldAcceleration = -constants.ElementaryCharge * m.Parameters.ConstEField / constants.ElectornMass

	m.sheathTimeOuterConstant = 2 * m.Parameters.CathodeFallLength * math.Sqrt(constants.ElectornMass/(2*constants.ElementaryCharge*m.Vc))
	m.sheathTimeInnerConstant = math.Sqrt(m.Parameters.CathodeFallPotential) / m.Parameters.CathodeFallLength
	m.constFieldTimeConstant = math.Sqrt(2*constants.ElectornMass/constants.ElementaryCharge) * (-m.Parameters.ConstEField)

	m.EFieldA = 2 * m.Parameters.CathodeFallPotential / (m.Parameters.CathodeFallLength * m.Parameters.CathodeFallLength)
	m.EFieldBb = -2*m.Parameters.CathodeFallPotential/m.Parameters.CathodeFallLength + m.Parameters.ConstEField
	m.VOffset = -m.Va + m.Parameters.CathodeFallLength*(0.5*m.EFieldA*m.Parameters.CathodeFallLength+m.EFieldBb)
	m.motionW = math.Sqrt(constants.ElectronChargeToMassRatio * m.EFieldA)
	m.motionK = -m.EFieldBb / m.EFieldA
}

func NewModel(parameters config.ModelParameters) *Model {
	m := Model{}
	m.Parameters = parameters
	m.InitVars()

	meanFreePath := 1. / (m.Parameters.CrossSectionsData().SurplusCrossSection() * m.Parameters.GasDensity)

	m.NumCells = 5 * int(m.Parameters.SimulationLength()/meanFreePath)
	m.NumCellsMu = int(2./parameters.MuDiscretizationStep) + 1
	m.NumCellsE = int((m.Va+m.Vc+5)/m.Parameters.EnergyDiscretizationStep) + 1

	m.XStep = m.Parameters.SimulationLength() / float64(m.NumCells)

	m.RateLookups = make(map[string][]float64)
	for p := range parameters.CrossSectionsData().Processes {
		if parameters.CrossSectionsData().Processes[p].Type == lxgata.IONIZATION {
			m.ionizationCS = &parameters.CrossSectionsData().Processes[p]
			m.RateLookups[string(lxgata.IONIZATION)] = make([]float64, m.NumCellsE+1)
			for e := range m.RateLookups[string(lxgata.IONIZATION)] {
				energy := (float64(e) + 0.5) * parameters.EnergyDiscretizationStep
				m.RateLookups[string(lxgata.IONIZATION)][e] = utils.EV2electronVelocity(energy) * parameters.CrossSectionsData().Processes[p].CrossSectionAt(energy) * parameters.GasDensity
			}
		} else {
			for substring := range parameters.RequireCollisionRelativeMargin {
				collisionName := string(parameters.CrossSectionsData().Processes[p].Type) + parameters.CrossSectionsData().Processes[p].Outcome
				if strings.Contains(collisionName, substring) {
					m.RateLookups[collisionName] = make([]float64, m.NumCellsE+1)
					for e := range m.RateLookups[string(lxgata.IONIZATION)] {
						energy := (float64(e) + 0.5) * parameters.EnergyDiscretizationStep
						m.RateLookups[collisionName][e] = utils.EV2electronVelocity(energy) * parameters.CrossSectionsData().Processes[p].CrossSectionAt(energy) * parameters.GasDensity
					}
				}
			}
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

	m.MeanFreePath = make([]utils.AggregationF, m.NumCells)

	numEnergyCells := int((m.Vc+m.Va+5+0.1)/m.Parameters.EnergyStep) + 10
	m.LookupEnergy = make([]float64, numEnergyCells*5)
	m.lookupMNumerator = make([]float64, numEnergyCells*5)
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
}

func (m *Model) nextCollisionTime(p *Particle, tupd chan<- TrajectoryUpdate) (collisionDescription *lxgata.Collision, throwOut, reversal bool, freePath utils.KahanSummable) {
	var totalCrossSectionPrimed float64
	totalCrossSectionPrimed = m.Parameters.CrossSectionsData().SurplusCrossSection()
	// if p.mu < 0 {
	// 	totalCrossSectionPrimed = m.Parameters.CrossSectionsData().MaxTotalCrossSectionOverRange(p.trajectory.radialEnergy, p.trajectory.totEnergy)
	// } else {
	// 	totalCrossSectionPrimed = m.Parameters.CrossSectionsData().MaxTotalCrossSectionOverRange(p.eKinetic, p.trajectory.totEnergy)
	// }
	maxCollisionFrequency := utils.EV2electronVelocity(p.trajectory.totEnergy) * totalCrossSectionPrimed * m.Parameters.GasDensity
	var nextCollisionTime, collisionEnergy float64
	lastSegmentPrograde := p.mu > 0
	collisionPosition := p.x
	var potentialAtCollisionPoint float64
	for collisionDescription == nil {
		nextCollisionTime = rand.ExpFloat64() / maxCollisionFrequency
		collisionPosition, lastSegmentPrograde = p.trajectory.PositionAfterTime(collisionPosition, nextCollisionTime, m, lastSegmentPrograde)
		if collisionPosition > m.Parameters.SimulationLength() {
			collisionPosition = m.Parameters.SimulationLength()
			collisionEnergy = p.trajectory.totEnergy + m.VfromL(collisionPosition)
			throwOut = true
			break
		}
		if collisionPosition < 0 {
			collisionPosition = 0
			collisionEnergy = p.trajectory.totEnergy + m.VfromL(collisionPosition)
			throwOut = true
			break
		}
		potentialAtCollisionPoint = m.VfromL(collisionPosition)
		collisionEnergy = p.trajectory.totEnergy + potentialAtCollisionPoint
		if collisionEnergy < p.trajectory.radialEnergy {
			collisionEnergy = p.trajectory.radialEnergy
		}
		collisionDescription = m.Parameters.CrossSectionsData().SampleWithNullCollision(collisionEnergy, totalCrossSectionPrimed)
	}
	if lastSegmentPrograde != (p.mu > 0) { //turnaround occured
		update1 := TrajectoryUpdate{
			startEnergy: p.eKinetic,
			endEnergy:   p.trajectory.radialEnergy,
			trajectory:  p.trajectory,
			origin:      p.origin,
		}
		reversal = true
		update2 := TrajectoryUpdate{
			startEnergy: p.trajectory.radialEnergy,
			endEnergy:   collisionEnergy,
			trajectory:  p.trajectory,
			origin:      p.origin,
		}
		// freePath.Add()
		tupd <- update1
		tupd <- update2
		p.mu = -p.mu
	} else {
		tupd <- TrajectoryUpdate{
			startEnergy: p.eKinetic,
			endEnergy:   collisionEnergy,
			trajectory:  p.trajectory,
			origin:      p.origin,
		}
	}
	if !throwOut {
		p.setEnergy(collisionEnergy, m, true, true)
	}
	return
}

type BoundaryConditionType int

const (
	None BoundaryConditionType = iota
	ArrivalAtCathode
	ArrivalAtSimEnd
	ArrivalAtTubeWall
)

func (m *Model) nextCollision(p *Particle, tupd chan<- TrajectoryUpdate) (collisionDescription *lxgata.Collision, reversal bool, throwOut BoundaryConditionType, freePath utils.KahanSummable) {
	R := -math.Log(1. - rand.Float64())

	if p.x < 0 {
		throwOut = ArrivalAtCathode
		return
	}
	if p.x > m.Parameters.SimulationLength() {
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
		if p.mu < 0 {
			if !alignedToEnergyGrid {
				nextCellIndex = int(p.eKinetic / m.Parameters.EnergyStep)
				lowEnergy, highEnergy = m.LookupEnergy[nextCellIndex], p.eKinetic

			} else {
				if currentCellIndex == 0 {
					lowEnergy, highEnergy = -m.Parameters.EnergyStep, 0
				} else {
					lowEnergy, highEnergy = m.LookupEnergy[currentCellIndex-1], m.LookupEnergy[currentCellIndex]
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
			if p.trajectory.totEnergy-lowEnergy > m.Vc+m.Va {
				arrivalAtCathode = true
				lowEnergy = p.trajectory.totEnergy - (m.Vc + m.Va)
			}
			segmentStartEnergy, segmentEndEnergy = highEnergy, lowEnergy
		} else {
			if !alignedToEnergyGrid {
				nextCellIndex = int(p.eKinetic/m.Parameters.EnergyStep) + 1
				lowEnergy, highEnergy = p.eKinetic, m.LookupEnergy[nextCellIndex]
			} else {
				// if currentCellIndex+1 >= len(m.LookupEnergy) {

				// }
				lowEnergy, highEnergy = m.LookupEnergy[currentCellIndex], m.LookupEnergy[currentCellIndex+1]
				nextCellIndex = currentCellIndex + 1
			}

			//check for gap end arrival
			if p.trajectory.totEnergy <= highEnergy {
				arrivalAtSimEnd = true
				highEnergy, highEnergyAligned = p.trajectory.totEnergy, false
			} else {
				highEnergyCellIndex, highEnergyAligned = nextCellIndex, true
			}
			segmentStartEnergy, segmentEndEnergy = lowEnergy, highEnergy
		}

		var M, G float64
		if highEnergyAligned {
			M = m.lookupMNumerator[highEnergyCellIndex] / -m.EFieldFromPotential(-(p.trajectory.totEnergy - highEnergy))
		} else {
			M = m.Parameters.GasDensity * m.Parameters.CrossSectionsData().TotalCrossSectionAt(highEnergy) * math.Sqrt(highEnergy) / -m.EFieldFromPotential(-(p.trajectory.totEnergy - highEnergy))
		}

		var higherVelocity, lowerVelocity float64
		if isVelocityCached {
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
		if G < R {
			R -= G
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
					delta = math.Sqrt(segmentStartEnergy-p.trajectory.radialEnergy) - R/(2.*M) // == sqrt(pColl - p.trajectory.eStar) i.e. abs(speed) equivalent along x axis
					if delta < 0 {
						R -= 2 * M * math.Sqrt(segmentStartEnergy-p.trajectory.radialEnergy)
						segmentEndEnergy = p.trajectory.radialEnergy
						reversalOccured = true
					}
					if math.IsNaN(delta) {
						panic("delta is NaN")
					}
				} else {
					delta = R/(2.*M) + math.Sqrt(segmentStartEnergy-p.trajectory.radialEnergy)
					if math.IsNaN(delta) {
						panic("delta is NaN")
					}
				}
				if !reversalOccured {
					collisionEnergy := math.FMA(delta, delta, p.trajectory.radialEnergy) // R = 2M[sqrt(p.e-p.trajectory.eStar) - sqrt(eColl - p.trajectory.eStar)]
					segmentEndEnergy = collisionEnergy
					p.setEnergy(collisionEnergy, m, true, true)
					if p.trajectory.totEnergy < collisionEnergy || m.Parameters.SimulationLength() < p.x {
						arrivalAtSimEnd = true
						segmentEndEnergy = p.trajectory.totEnergy
					} else if p.x < 0 {
						arrivalAtCathode = true
						segmentEndEnergy = p.trajectory.totEnergy - (m.Vc + m.Va)
					} else {
						if !m.Parameters.Volumetric || p.y*p.y+p.z*p.z < m.tubeRadius2 {
							var totalCrossSectionPrimed = M * math.Abs(m.EFieldFromL(p.x)) / (math.Sqrt(collisionEnergy) * m.Parameters.GasDensity)
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
			p.x = 0
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
				p.setEnergy(p.trajectory.totEnergy, m, true, false)
				alignedToEnergyGrid = false
			} else {
				if m.Parameters.AnodeBackscatteringCoefficient0 != 0 {
					backscatteringCoefficient := m.Parameters.AnodeBackscatteringCoefficient0 * math.Exp(m.Parameters.AnodeBackscatteringCoefficientB*(1-p.mu))

					p.setEnergy(p.trajectory.totEnergy, m, true, true)
					if rand.Float64() < backscatteringCoefficient && p.y*p.y+p.z*p.z < m.cathodeRadius2 {
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
					p.x = m.Parameters.SimulationLength()
					p.eKinetic = p.trajectory.totEnergy
					p.mu = math.Sqrt(p.getAxialEnergy() / p.trajectory.totEnergy)
					break
				}

			}
		} else if reversalOccured {
			p.mu = +0.
			p.eKinetic = p.trajectory.radialEnergy
			alignedToEnergyGrid = false
			// pathStartEnergy = p.trajectory.radialEnergy
			reversal = true
		} else {
			alignedToEnergyGrid = true
		}
		if m.Parameters.Volumetric && !math.IsInf(m.tubeRadius2, 0) { //&& currentCellIndex == wallCollisionEnergyStep
			throwOut = ArrivalAtTubeWall
			break
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
			}
			tupd <- TrajectoryUpdate{
				startEnergy: p.trajectory.radialEnergy,
				endEnergy:   segmentEndEnergy,
				trajectory:  p.trajectory,
				origin:      p.origin,
			}
		} else {
			tupd <- TrajectoryUpdate{
				startEnergy: pathStartEnergy,
				endEnergy:   segmentEndEnergy,
				trajectory:  p.trajectory,
				origin:      p.origin,
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
}

type MeanFreePathUpdate struct {
	freePath                     utils.KahanSummable
	fillStartIndex, fillEndIndex int
}

func (m *Model) Run(electronsToSimulate func(*Model) int) {
	CollisionAtCellPerAvalanche := make(map[lxgata.CollisionType][][]int)
	for process := range m.CollisionAtCell {
		CollisionAtCellPerAvalanche[process] = make([][]int, m.NumCells+1)
	}
	DetailedCollisionAtCell := make(map[string][][]int)
	for process := range m.DetailedCollisionAtCell {
		DetailedCollisionAtCell[process] = make([][]int, m.NumCells+1)
	}

	TimeAtCellPerAvalanche := make([][]utils.KahanSummable, m.NumCells+1)

	var energyDeposition float64
	var numberOfEDEvents int

	nElectrons := electronsToSimulate(m)
	for nElectrons > 0 {
		m.PurgeMetrics()

		computeFlow := make(chan *Particle, nElectrons*1000)
		var computeWg, stateWg sync.WaitGroup

		collFlow := make(chan CollisionEvent, 1000*m.Parameters.NElectrons)
		stateWg.Add(1)
		go func() {
			for collision := range collFlow {
				if collision.x < m.NumCells {
					CollisionAtCellPerAvalanche[collision.collType][collision.x][collision.origin-m.TotalElectronsEmittedOnCathode]++
					collName := string(collision.collType) + collision.outcome
					if _, exist := DetailedCollisionAtCell[collName]; exist {
						DetailedCollisionAtCell[collName][collision.x][collision.origin-m.TotalElectronsEmittedOnCathode]++
					}
				}
			}
			stateWg.Done()
		}()

		var trajFlow chan TrajectoryUpdate
		if m.Parameters.CellTimeWeighting || m.Parameters.CalculateDistribution {
			for x := range TimeAtCellPerAvalanche {
				TimeAtCellPerAvalanche[x] = make([]utils.KahanSummable, nElectrons)
			}

			trajFlow = make(chan TrajectoryUpdate, 100000)
			stateWg.Add(1)
			go func() {
				if m.Parameters.CellTimeWeighting || m.Parameters.CalculateDistribution {
					for tUpdate := range trajFlow {

						xStart := m.LfromV(tUpdate.startEnergy - tUpdate.trajectory.totEnergy)
						xEnd := m.LfromV(tUpdate.endEnergy - tUpdate.trajectory.totEnergy)
						xStart, xEnd = min(xStart, xEnd), max(xStart, xEnd)
						xStartIndex := int(math.Ceil(xStart / m.XStep))
						xEndIndex := int(math.Floor(xEnd / m.XStep))
						for xIndex := xStartIndex; xIndex <= xEndIndex; xIndex++ {
							energy := tUpdate.trajectory.totEnergy + m.lookUpPotentialAtCellNode[xIndex]
							eIndex := int(energy / m.Parameters.EnergyDiscretizationStep)
							mu := math.Copysign(math.Sqrt((energy-tUpdate.trajectory.radialEnergy)/energy), tUpdate.endEnergy-tUpdate.startEnergy)
							muIndex := int((mu + 1 + 0.5*m.Parameters.MuDiscretizationStep) / m.Parameters.MuDiscretizationStep)
							if xIndex < m.NumCells && eIndex < m.NumCellsE && muIndex < m.NumCellsMu {
								m.MeanEnergy[xIndex].AppendF(energy)
								m.DistributionXEMu[xIndex][eIndex][muIndex]++
							}
						}
					}
				}
				stateWg.Done()
			}()
		}

		for _, process := range m.Parameters.CrossSectionsData().GetTypes() {
			for x := range CollisionAtCellPerAvalanche[process] {
				CollisionAtCellPerAvalanche[process][x] = make([]int, nElectrons)
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

		electronsReturned := make([]int, nElectrons)

		for origin := range nElectrons {
			particle := m.newParticle(origin + m.TotalElectronsEmittedOnCathode)
			if m.Parameters.CalculateDistribution {
				m.DistributionXEMu[0][int(particle.eKinetic/m.Parameters.EnergyDiscretizationStep)][int((particle.mu+1)/m.Parameters.MuDiscretizationStep)]++
			}
			computeWg.Add(1)
			computeFlow <- &particle
		}

		status := []string{"//", "==", "\\\\", "||"}
		for range m.Parameters.Threads() {
			go func() {
				counter := 0
				for particlePtr := range computeFlow {
					if !m.Parameters.SupressSpinner {
						counter++
						print("\r" + status[counter&0b11])
					}

					lowerEnergyThreshold := m.Parameters.LowerEnergyThreshold()
					for lowerEnergyThreshold < particlePtr.trajectory.totEnergy && (!m.Parameters.ReturningElectrons || particlePtr.trajectory.totEnergy+m.VfromL(0)-particlePtr.trajectory.radialEnergy > 0) {
						fillStart := particlePtr.x
						collision, reversal, throwOut, freePath := m.nextCollision(particlePtr, trajFlow /*, stateflow*/)
						if m.Parameters.MeanFreePath {
							if reversal {
								fillStart = particlePtr.trajectory.getTurnaroundX(m)
							}
							fillEnd := particlePtr.x
							fillStart, fillEnd = min(fillStart, fillEnd), max(fillStart, fillEnd)
							freePathFlow <- MeanFreePathUpdate{
								freePath:       freePath,
								fillStartIndex: int(math.Ceil(fillStart / m.XStep)),
								fillEndIndex:   int(math.Floor(fillEnd / m.XStep)),
							}
						}

						if throwOut != None {
							if throwOut == ArrivalAtCathode {
								electronsReturned[particlePtr.origin]++
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
								ejected := *particlePtr
								ejected.ejectedFromIonization = true

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

								ejected.redirect(cosChiEjected, math.Cos(phi+math.Pi), m)
								if ejected.trajectory.totEnergy >= lowerEnergyThreshold {
									computeWg.Add(1)
									computeFlow <- &ejected
								}

							case lxgata.ATTACHMENT:
							case lxgata.EFFECTIVE:
							case lxgata.EXCITATION: //energy is lost with threshold decrement
							case lxgata.ROTATION:
							case lxgata.DEEXCITATION: //energy is gained with (negative) threshold decrement

							default:
								panic(fmt.Sprintf("unexpected lxgata.CollisionType: %#v", collision.Type))
							}
							// if !m.Parameters.Volumetric || (particlePtr.y*particlePtr.y+particlePtr.z*particlePtr.z < m.cathodeRadius2) || particlePtr.x < m.Parameters.CathodeFallLength {

							collFlow <- CollisionEvent{
								x:          int(particlePtr.x / m.XStep),
								energyLoss: energyLoss,
								collType:   collision.Type,
								origin:     particlePtr.origin,
								outcome:    collision.Outcome,
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
								break
							}
							particlePtr.redirect(cosChiScattered, math.Cos(phi), m)
						}
					}
					computeWg.Done()
				}
			}()
		}
		computeWg.Wait()

		close(computeFlow)
		close(collFlow)
		if m.Parameters.CellTimeWeighting || m.Parameters.CalculateDistribution {
			close(trajFlow)
		}

		if m.Parameters.MeanFreePath {
			close(freePathFlow)
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

		PerAvalancheCollisionSumsUpToCell := make([]int, nElectrons)
		var PerAvalancheTimeSumsUpToCell []utils.KahanSummable
		if m.Parameters.CellTimeWeighting {
			PerAvalancheTimeSumsUpToCell = make([]utils.KahanSummable, nElectrons)
		}
		for x := range m.IonizationsSumUpToCell {
			for electron := range nElectrons {
				PerAvalancheCollisionSumsUpToCell[electron] += CollisionAtCellPerAvalanche[lxgata.IONIZATION][x][electron]
				if m.Parameters.CellTimeWeighting {
					PerAvalancheTimeSumsUpToCell[electron].Add(TimeAtCellPerAvalanche[x][electron].Val())
				}
			}
			m.IonizationsSumUpToCell[x].Update(PerAvalancheCollisionSumsUpToCell)

		}
		m.ElectronsReturned.Update(electronsReturned)
		m.TotalElectronsEmittedOnCathode += nElectrons
		nElectrons = electronsToSimulate(m)
	}

	if numberOfEDEvents > 0 {
		m.MeanElectronEnergyAtAnode = energyDeposition / float64(numberOfEDEvents)
	}
	print("\r")
}
