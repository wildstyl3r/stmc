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
	inverseConstEField       float64
	inverseNegativeVc        float64
	inverseCathodeFallLength float64
	sheathTimeOuterConstant  float64
	sheathTimeInnerConstant  float64
	constFieldTimeConstant   float64
	constFieldAcceleration   float64
	timeMotionK              float64
	timeMotionW              float64
	timeMotionP              float64

	electronEnergyInEvToVelocityFactor float64
	// boundaryThickness        float64
	cathodeRadius2, tubeRadius2 float64
	// boundaryRadius2          float64

	Va float64 // additional voltage to avoid numeric negative energy beyond cathode fall

	NumCells, NumCellsMu, NumCellsE int

	XStep float64

	lookUpPotentialAtCellCenter []float64
	LookupEnergy                []float64
	// LookupInverseVelocity     []float64
	// GridInverseMu           []float64
	lookupMNumerator []float64

	DistributionXEMu [][][]float64
	DistributionXE   [][]float64
	// DistributionXE   [][][]float64

	// TimeAtCell      [][]utils.AggregationF //(x, energy) -- average over avalanches (on average update null times must be included)
	RateLookups  map[string][]float64
	ionizationCS *lxgata.Collision
	// timeAtCell      []utils.KahanSummable
	TimeAtCell      []utils.AggregationF
	CollisionAtCell map[lxgata.CollisionType][]utils.Aggregation ////utils.MutableMap[lxgata.CollisionType, utils.AggregatedDistribution]
	// CollisionOutside        utils.MutableMap[lxgata.CollisionType, utils.AggregatedDistribution]
	DetailedCollisionAtCell map[string][]utils.Aggregation //utils.MutableMap[string, utils.AggregatedDistribution]
	// MeanEnergy              utils.AggregatedDistribution
	IonizationsSumUpToCell []utils.Aggregation //utils.AggregatedDistribution
	// WallLossAtCell         utils.AggregatedDistribution
	// EnergyLossByProcess map[lxgata.CollisionType][]float64

	// lookupCosinesAtAngleCellBounds  []float64
	// LookupCosinesAtAngleCellCenters []float64

	TotalElectronsEmittedOnCathode int

	MeanElectronEnergyAtAnode float64

	OutOfEnergyAtCell []int
	DataHub           *DataHubType
}

func (m *Model) CalculateStaticTimeMotionConstants() {
	m.timeMotionK = m.Parameters.ConstEField*m.Parameters.CathodeFallLength + m.Parameters.CathodeFallPotential
	m.timeMotionW = math.Sqrt(2*constants.ElementaryCharge*m.timeMotionK/constants.ElectornMass) / m.Parameters.CathodeFallLength
	m.timeMotionP = 0.5 * (1. + m.Parameters.CathodeFallPotential/m.timeMotionK) * m.Parameters.CathodeFallLength
}

func NewModel(parameters config.ModelParameters) *Model {
	m := Model{}
	m.Parameters = parameters
	m.Vc = math.Abs(m.Parameters.CathodeFallPotential)
	m.Va = math.Abs(m.Parameters.ConstEField * (m.Parameters.SimulationLength() - m.Parameters.CathodeFallLength))
	m.inverseConstEField = 1. / m.Parameters.ConstEField
	m.inverseNegativeVc = -1. / m.Vc
	m.inverseCathodeFallLength = 1. / m.Parameters.CathodeFallLength
	m.cathodeRadius2 = parameters.CathodeRadius * parameters.CathodeRadius
	m.tubeRadius2 = parameters.TubeRadius * parameters.TubeRadius
	m.constFieldAcceleration = -constants.ElementaryCharge * parameters.ConstEField / constants.ElectornMass

	m.sheathTimeOuterConstant = 2 * parameters.CathodeFallLength * math.Sqrt(constants.ElectornMass/(2*constants.ElementaryCharge*m.Vc))
	m.sheathTimeInnerConstant = math.Sqrt(parameters.CathodeFallPotential) / parameters.CathodeFallLength
	m.constFieldTimeConstant = math.Sqrt(2*constants.ElectornMass/constants.ElementaryCharge) * (-parameters.ConstEField)

	m.electronEnergyInEvToVelocityFactor = math.Sqrt(2 * constants.ElectronChargeToMassRatio)

	m.CalculateStaticTimeMotionConstants()

	meanFreePath := 1. / (m.Parameters.CrossSectionsData().SurplusCrossSection() * m.Parameters.GasDensity)
	// if m.Parameters.Verbose() {
	// 	fmt.Printf("Mean free path: %f\n", meanFreePath)
	// }

	m.NumCells = 5 * int(m.Parameters.SimulationLength()/meanFreePath)
	m.NumCellsMu = int(2./parameters.MuDiscretizationStep) + 1
	m.NumCellsE = int((m.Va+m.Vc+5)/m.Parameters.EnergyDiscretizationStep) + 1
	// if m.Parameters.Verbose() {
	// 	fmt.Printf("xCells: %d;\n", m.NumCells)
	// }

	m.XStep = m.Parameters.SimulationLength() / float64(m.NumCells)
	// m.MuStep = m.Parameters.MuStep    // eV
	// m.NumMuCells = int(2. / m.MuStep)

	// radianAngleStep := parameters.AngleStep / 360. * 2. * math.Pi
	// numAngleCells := int(360. / parameters.AngleStep)
	// m.lookupCosinesAtAngleCellBounds, m.LookupCosinesAtAngleCellCenters = make([]float64, numAngleCells), make([]float64, numAngleCells)
	// for a := range numAngleCells {
	// 	angle := float64(a) * radianAngleStep
	// 	nextAngle := float64(a+1) * radianAngleStep
	// 	m.lookupCosinesAtAngleCellBounds[a] = math.Cos(angle)
	// 	m.LookupCosinesAtAngleCellCenters[a] = math.Cos(0.5 * (angle + nextAngle))
	// }

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
	// m.TimeAtCell = make([][]utils.AggregationF, m.NumCells+1)
	// for x := range m.TimeAtCell {
	// 	m.TimeAtCell[x] = make([]utils.AggregationF, m.NumCellsE+1)
	// }

	m.CollisionAtCell = make(map[lxgata.CollisionType][]utils.Aggregation) //utils.NewMutableMap[lxgata.CollisionType, utils.AggregatedDistribution]() //make(map[lxgata.CollisionType]utils.AggregatedDistribution)
	// m.CollisionOutside = utils.NewMutableMap[lxgata.CollisionType, utils.AggregatedDistribution]() //make(map[lxgata.CollisionType]utils.AggregatedDistribution)
	// m.DetailedCollisionAtCell = make(map[string][]utils.AggregationF) //utils.NewMutableMap[string, utils.AggregatedDistribution]() //make(map[string]utils.AggregatedDistribution)
	// m.EnergyLossByProcess = make(map[lxgata.CollisionType][]float64)
	m.IonizationsSumUpToCell = make([]utils.Aggregation, m.NumCells+1) //utils.NewAggregatedDistribution(m.NumCells + 1)
	// m.MeanEnergy = utils.NewAggregatedDistribution(m.NumCells + 1)
	for _, process := range parameters.CrossSectionsData().GetTypes() {
		// m.EnergyLossByProcess[process] = make([]float64, m.NumCells+1)
		m.CollisionAtCell[process] = make([]utils.Aggregation, m.NumCells+1) //utils.NewAggregatedDistribution(m.NumCells + 1)
		// 	*m.CollisionOutside.GetPointer(process) = utils.NewAggregatedDistribution(m.NumCells + 1)
	}
	m.DetailedCollisionAtCell = make(map[string][]utils.Aggregation) //utils.NewMutableMap[string, utils.AggregatedDistribution]() //make(map[string]utils.AggregatedDistribution)
	for p := range parameters.CrossSectionsData().Processes {
		for substring := range parameters.RequireCollisionRelativeMargin {
			collisionName := string(parameters.CrossSectionsData().Processes[p].Type) + parameters.CrossSectionsData().Processes[p].Outcome
			if strings.Contains(collisionName, substring) {
				// *m.DetailedCollisionAtCell.GetPointer(collisionName) = utils.NewAggregatedDistribution(m.NumCells + 1)
				m.DetailedCollisionAtCell[collisionName] = make([]utils.Aggregation, m.NumCells+1)
			}
		}
	}
	// m.WallLossAtCell = utils.NewAggregatedDistribution(m.NumCells + 1)

	m.lookUpPotentialAtCellCenter = make([]float64, m.NumCells+2)
	for cellNode := range m.lookUpPotentialAtCellCenter {
		m.lookUpPotentialAtCellCenter[cellNode] = m.VfromL((float64(cellNode) + 0.5) * m.XStep)
	}

	m.OutOfEnergyAtCell = make([]int, m.NumCells)

	numEnergyCells := int((m.Vc+m.Va+5+0.1)/m.Parameters.EnergyStep) + 1
	m.LookupEnergy = make([]float64, numEnergyCells)
	m.lookupMNumerator = make([]float64, numEnergyCells)
	for gridNode := range m.LookupEnergy {
		m.LookupEnergy[gridNode] = float64(gridNode) * m.Parameters.EnergyStep
		m.lookupMNumerator[gridNode] = m.Parameters.GasDensity * m.Parameters.CrossSectionsData().TotalCrossSectionAt(m.LookupEnergy[gridNode]) * math.Sqrt(m.LookupEnergy[gridNode])

	}
	// m.GridInverseMu = make([]float64, int(1./m.MuStep)+1)
	// for muNode := range m.GridInverseMu {
	// 	m.GridInverseMu[muNode] = 1. / (m.MuStep * (float64(muNode) + 0.5))
	// }

	if m.Parameters.CalculateDistribution {
		m.DistributionXEMu = make([][][]float64, m.NumCells)
		for x := range m.DistributionXEMu {
			m.DistributionXEMu[x] = make([][]float64, m.NumCellsE)
			for e := range m.DistributionXEMu[x] {
				m.DistributionXEMu[x][e] = make([]float64, m.NumCellsMu)
			}
		}
		m.DistributionXE = make([][]float64, m.NumCells)
		for x := range m.DistributionXE {
			m.DistributionXE[x] = make([]float64, m.NumCellsE)
		}
		// m.DistributionXE = make([][][]float64, m.NumCells)
		// for x := range m.DistributionXE {
		// 	m.DistributionXE[x] = make([][]float64, m.NumCellsE)
		// 	for e := range m.DistributionXE[x] {
		// 		m.DistributionXE[x][e] = make([]float64, 0)
		// 	}
		// }
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

func (m *Model) ConstEFieldTime(axialEnergy1_eV, axialEnergy2_eV float64) float64 { //time to get from x1 to x2 in the constant field region
	return m.constFieldTimeConstant * math.Abs(math.Sqrt(axialEnergy2_eV)-math.Sqrt(axialEnergy1_eV))
}

func (m *Model) collisionSelector(eKinetic, x, M float64) *lxgata.Collision {
	var totalCrossSectionPrimed = M * math.Abs(m.EFieldFromL(x)) / math.Sqrt(eKinetic) / m.Parameters.GasDensity
	return m.Parameters.CrossSectionsData().SampleWithNullCollision(eKinetic, totalCrossSectionPrimed)
}

type TrajectoryUpdate struct {
	startEnergy, endEnergy float64
	trajectory             TrajectoryConstants
	origin                 int
}

func (m *Model) timeToWall(p *Particle) float64 {
	if math.IsInf(m.Parameters.TubeRadius, 0) {
		return math.Inf(-1)
	}
	vR := utils.EV2electronVelocity(p.trajectory.radialEnergy)
	normY, normZ, normR := p.y/vR, p.z/vR, m.Parameters.TubeRadius/vR
	halfB := normY*p.sinEta + normZ*p.cosEta
	D := halfB*halfB - (normY*normY + normZ*normZ - normR*normR)
	if D >= 0 {
		s1 := -halfB - math.Sqrt(D)
		s2 := -halfB + math.Sqrt(D)
		if s1 > 0 {
			return s1
		} else if s2 > 0 {
			return s2
		} else {
			return math.Inf(-1)
		}
	} else {
		return math.Inf(1)
	}
}

// func (m *Model) nextCollisionTime(p *Particle, tupd chan<- TrajectoryUpdate) (collisionDescription *lxgata.Collision, throwOut bool) {
// 	var totalCrossSectionPrimed float64
// 	if p.mu < 0 {
// 		totalCrossSectionPrimed = m.Parameters.CrossSectionsData().MaxTotalCrossSectionOverRange(p.trajectory.radialEnergy, p.trajectory.totEnergy)
// 	} else {
// 		totalCrossSectionPrimed = m.Parameters.CrossSectionsData().MaxTotalCrossSectionOverRange(p.eKinetic, p.trajectory.totEnergy)
// 	}
// 	maxCollisionFrequency := utils.EV2electronVelocity(p.trajectory.totEnergy) * totalCrossSectionPrimed * m.Parameters.GasDensity
// 	var nextCollisionTime, collisionEnergy, collisionPosition float64
// 	for collisionDescription == nil {
// 		nextCollisionTime = rand.ExpFloat64() / maxCollisionFrequency
// 		collisionPosition = localTimeMotionConstants.GetPositionAndEnergy(nextCollisionTime)
// 		collisionEnergy = p.trajectory.totEnergy + m.VfromL(collisionPosition)
// 		collisionDescription = m.Parameters.CrossSectionsData().SampleWithNullCollision(collisionEnergy, totalCrossSectionPrimed)
// 	}
// 	tupd <- TrajectoryUpdate{
// 		startEnergy:  p.eKinetic,
// 		endEnergy:    collisionEnergy,
// 		totalEnergy:  p.trajectory.totEnergy,
// 		radialEnergy: p.trajectory.eStar,
// 		origin:       p.origin,
// 	}
// 	p.setEnergy(collisionEnergy, m, true, true)
// }

// DO NOT USE
func (m *Model) nextCollisionSpatial(p *Particle, tupd chan<- TrajectoryUpdate) (collisionDescription *lxgata.Collision, throwOut bool) {
	pathStartingPoint := p.x
	var collisionEnergy float64
	for collisionDescription == nil {
		var totalCrossSectionPrimed float64
		if p.mu < 0 {
			totalCrossSectionPrimed = m.Parameters.CrossSectionsData().MaxTotalCrossSectionOverRange(p.trajectory.radialEnergy, p.trajectory.totEnergy)
		} else {
			totalCrossSectionPrimed = m.Parameters.CrossSectionsData().MaxTotalCrossSectionOverRange(p.eKinetic, p.trajectory.totEnergy)
		}
		R := -math.Log(1. - rand.Float64())
		// var reversalOccured, arrivalAtCathode, arrivalAtSimEnd, arrivalAtTubeWall bool
		pathLength := R / (totalCrossSectionPrimed * m.Parameters.GasDensity)
		var collisionPosition float64
		if p.mu > 0 {
			collisionPosition = p.x + pathLength
		} else {
			turnaroundPoint := p.trajectory.getTurnaroundX(m)
			if p.x-turnaroundPoint > pathLength {
				collisionPosition = p.x - pathLength
			} else {
				if turnaroundPoint < 0 {
					throwOut = true
					tupd <- TrajectoryUpdate{
						startEnergy: p.trajectory.totEnergy + m.VfromL(pathStartingPoint),
						endEnergy:   p.trajectory.radialEnergy,
						trajectory:  p.trajectory,
						origin:      p.origin,
					}
					break
				}
				tupd <- TrajectoryUpdate{
					startEnergy: p.trajectory.totEnergy + m.VfromL(pathStartingPoint),
					endEnergy:   p.trajectory.radialEnergy,
					trajectory:  p.trajectory,
					origin:      p.origin,
				}
				pathStartingPoint = turnaroundPoint
				p.mu = -p.mu
				collisionPosition = turnaroundPoint + (pathLength - (p.x - turnaroundPoint))
			}
		}
		if collisionPosition > m.Parameters.SimulationLength() {
			// arrivalAtSimEnd = true
			p.x = m.Parameters.SimulationLength()
			p.eKinetic = p.trajectory.totEnergy
			p.mu = math.Sqrt(p.getAxialEnergy() / p.trajectory.totEnergy)
			throwOut = true
			break
		} else if collisionPosition < 0 {
			// arrivalAtCathode = true
			throwOut = true
			break
		} else {
			collisionEnergy = p.trajectory.totEnergy + m.VfromL(collisionPosition)
			collisionDescription = m.Parameters.CrossSectionsData().SampleWithNullCollision(collisionEnergy, totalCrossSectionPrimed)
			tupd <- TrajectoryUpdate{
				startEnergy: p.trajectory.totEnergy + m.VfromL(pathStartingPoint),
				endEnergy:   collisionEnergy,
				trajectory:  p.trajectory,
				origin:      p.origin,
			}
			p.setEnergy(collisionEnergy, m, true, true)
		}
	}
	return
}

func (m *Model) nextCollision(p *Particle, tupd chan<- TrajectoryUpdate) (collisionDescription *lxgata.Collision, throwOut bool) {
	// wallCollisionEnergyStep := -1
	// var tw, xw, ew float64
	// if m.Parameters.Volumetric {
	// 	tw, xw = m.WallIntersection(p)
	// 	if !math.IsInf(tw, 0) {
	// 		ew = p.trajectory.totEnergy + m.VfromL(xw)
	// 		wallCollisionEnergyStep = int(ew / m.Parameters.EnergyStep)
	// 	}
	// }

	R := utils.NewKahanSummable(-math.Log(1. - rand.Float64()))

	if !(0 <= p.x && p.x < m.Parameters.SimulationLength()) {
		throwOut = true
		return
	}

	var currentCellIndex int
	var cachedVelocity float64
	var isVelocityCached = false

	alignedToEnergyGrid := false

	pathStartEnergy := p.eKinetic

	for {
		var lowEnergy, segmentStartEnergy, segmentEndEnergy, highEnergy float64
		var nextCellIndex, highEnergyCellIndex int
		var reversalOccured, arrivalAtCathode, arrivalAtSimEnd, arrivalAtTubeWall, highEnergyAligned bool
		if p.mu < 0 {
			if !alignedToEnergyGrid {
				nextCellIndex = int(p.eKinetic / m.Parameters.EnergyStep)
				lowEnergy, highEnergy = m.LookupEnergy[nextCellIndex], p.eKinetic

			} else {
				lowEnergy, highEnergy = m.LookupEnergy[currentCellIndex-1], m.LookupEnergy[currentCellIndex]
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
		if G < R.Val() {
			R.Add(-G)
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
					p.setEnergy(collisionEnergy, m, true, true)
					if p.trajectory.totEnergy < collisionEnergy || m.Parameters.SimulationLength() < p.x {
						arrivalAtSimEnd = true
						segmentEndEnergy = p.trajectory.totEnergy
					} else if p.x < 0 {
						arrivalAtCathode = true
						segmentEndEnergy = p.trajectory.totEnergy - (m.Vc + m.Va)
					} else {
						if !m.Parameters.Volumetric || p.y*p.y+p.z*p.z < m.tubeRadius2 {
							collisionDescription = m.collisionSelector(collisionEnergy, p.x, M)
							collisionOccured = true
						} else {
							arrivalAtTubeWall = true
						}
					}
				}
			}
		}

		if collisionOccured || reversalOccured || arrivalAtCathode || arrivalAtSimEnd || arrivalAtTubeWall {
			tupd <- TrajectoryUpdate{
				startEnergy: pathStartEnergy,
				endEnergy:   segmentEndEnergy,
				trajectory:  p.trajectory,
				origin:      p.origin,
			}
			// m.getFlowBetweenEnergies(pathStartEnergy, segmentEndEnergy, p.trajectory.totEnergy, p.trajectory.eStar, p.origin, tupd)
		}

		if collisionOccured {
			return
		}

		if arrivalAtCathode {
			p.x = 0
			throwOut = true
			return
		}

		// if arrivalAtTubeWall {
		// 	p.x = xw
		// 	var tubeBackscatter0 float64 //tubeBackscatterB
		// 	if m.Parameters.ParallelPlaneHollowCathode {
		// 		tubeBackscatter0 = m.Parameters.AnodeBackscatteringCoefficient0
		// 		// tubeBackscatterB = m.Parameters.AnodeBackscatteringCoefficientB
		// 	} else {
		// 		tubeBackscatter0 = m.Parameters.TubeBackscatteringCoefficient0
		// 		// tubeBackscatterB = m.Parameters.TubeBackscatteringCoefficientB
		// 	}
		// 	if tubeBackscatter0 != 0 {
		// 		// backscatteringCoefficient := tubeBackscatter0 * math.Exp(tubeBackscatterB*(1-math.Cos(incid)))
		// 	} else {
		// 		throwOut = true
		// 		return
		// 	}
		// }

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
						throwOut = true
						return
					}
				} else {
					throwOut = true
					p.x = m.Parameters.SimulationLength()
					p.eKinetic = p.trajectory.totEnergy
					p.mu = math.Sqrt(p.getAxialEnergy() / p.trajectory.totEnergy)
					return
				}

			}
		} else if reversalOccured {
			p.mu = +0.
			p.eKinetic = p.trajectory.radialEnergy
			alignedToEnergyGrid = false
			pathStartEnergy = p.trajectory.radialEnergy
		} else {
			alignedToEnergyGrid = true
		}
		if m.Parameters.Volumetric && !math.IsInf(m.tubeRadius2, 0) { //&& currentCellIndex == wallCollisionEnergyStep
			throwOut = true
			return
		}
		currentCellIndex = nextCellIndex
	}
}

type CollisionEvent struct {
	x          int
	energyLoss float64
	collType   lxgata.CollisionType
	outcome    string
	origin     int
}

type WallLossEvent struct {
	x      int
	absVx  float64
	origin int
}

func (m *Model) Run(electronsToSimulate func(*Model) int) {
	CollisionAtCellPerAvalanche := make(map[lxgata.CollisionType][][]int)
	// CollisionOutside_perElectron := make(map[lxgata.CollisionType][]utils.WeightedDistribution)
	for process := range m.CollisionAtCell {
		CollisionAtCellPerAvalanche[process] = make([][]int, m.NumCells+1)
		// CollisionOutside[process] = make([][]int, m.NumCells+1)
	}
	DetailedCollisionAtCell := make(map[string][][]int)
	for process := range m.DetailedCollisionAtCell {
		DetailedCollisionAtCell[process] = make([][]int, m.NumCells+1)
	}
	// var WallLossAtCell []utils.WeightedDistribution

	TimeAtCellPerAvalanche := make([][]utils.KahanSummable, m.NumCells+1)

	var energyDeposition float64
	var numberOfEDEvents int

	nElectrons := electronsToSimulate(m)
	for nElectrons > 0 {
		m.PurgeMetrics()

		computeFlow := make(chan *Particle, nElectrons*100)
		var computeWg, stateWg sync.WaitGroup

		energyDepositionFlow := make(chan float64, 10*m.Parameters.NElectrons)
		stateWg.Add(1)
		go func() {
			for deposition := range energyDepositionFlow {
				energyDeposition += deposition
				numberOfEDEvents += 1
			}
			stateWg.Done()
		}()

		collFlow := make(chan CollisionEvent, 10*m.Parameters.NElectrons)
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

		for x := range TimeAtCellPerAvalanche {
			TimeAtCellPerAvalanche[x] = make([]utils.KahanSummable, nElectrons)
		}

		trajFlow := make(chan TrajectoryUpdate, 1000)
		stateWg.Add(1)
		go func() {
			for tUpdate := range trajFlow {
				if m.Parameters.CellTimeWeighting {
					energyMin, energyMax := min(tUpdate.startEnergy, tUpdate.endEnergy), max(tUpdate.startEnergy, tUpdate.endEnergy)
					xMin, xMax := max(0, m.LfromV(-(tUpdate.trajectory.totEnergy-energyMin))), m.LfromV(-(tUpdate.trajectory.totEnergy - energyMax))
					xMinCrossingIndex, xMaxCrossingIndex := int(math.Ceil(xMin/m.XStep)), int(math.Floor(xMax/m.XStep))
					if xMinCrossingIndex > xMaxCrossingIndex {
						TimeAtCellPerAvalanche[xMaxCrossingIndex][tUpdate.origin-m.TotalElectronsEmittedOnCathode].Add(tUpdate.trajectory.TimeBetweenPositionsNoReversal(xMin, xMax, m))
					} else {
						if xMinCrossingIndex < 0 {
							panic("i<0")
						}
						if xMinCrossingIndex != 0 {
							TimeAtCellPerAvalanche[xMinCrossingIndex-1][tUpdate.origin-m.TotalElectronsEmittedOnCathode].Add(tUpdate.trajectory.TimeBetweenPositionsNoReversal(xMin, float64(xMinCrossingIndex)*m.XStep, m))
						}
						{
							for x := xMinCrossingIndex; x < xMaxCrossingIndex; x++ {
								TimeAtCellPerAvalanche[x][tUpdate.origin-m.TotalElectronsEmittedOnCathode].Add(tUpdate.trajectory.TimeBetweenPositionsNoReversal(float64(x)*m.XStep, float64(x+1)*m.XStep, m))
								// m.timeAtCell[x].Add(tUpdate.trajectory.TimeBetweenPositionsNoReversal(float64(x)*m.XStep, float64(x+1)*m.XStep, m))
							}
						}
						if xMaxCrossingIndex != len(TimeAtCellPerAvalanche) {
							TimeAtCellPerAvalanche[xMaxCrossingIndex][tUpdate.origin-m.TotalElectronsEmittedOnCathode].Add(tUpdate.trajectory.TimeBetweenPositionsNoReversal(float64(xMaxCrossingIndex)*m.XStep, xMax, m))
						}
					}
				}
			}
			// for tUpdate := range trajFlow {
			// 	energyMin, energyMax := min(tUpdate.startEnergy, tUpdate.endEnergy), max(tUpdate.startEnergy, tUpdate.endEnergy)
			// 	xMin, xMax := m.LfromV(-(tUpdate.totalEnergy - energyMin)), m.LfromV(-(tUpdate.totalEnergy - energyMax))
			// 	xMinCrossingIndex, xMaxCrossingIndex := int(math.Ceil(xMin/m.XStep)), int(math.Floor(xMax/m.XStep))
			// 	if xMinCrossingIndex < 0 {
			// 		panic("i<0")
			// 	}
			// 	xList := make([]float64, xMaxCrossingIndex-xMinCrossingIndex+1)
			// 	{
			// 		for x := xMinCrossingIndex; x <= xMaxCrossingIndex; x++ {
			// 			xList[x-xMinCrossingIndex] = float64(x) * m.XStep
			// 		}
			// 	}
			// 	eListFromX := make([]float64, len(xList))
			// 	for x := range xList {
			// 		eListFromX[x] = tUpdate.totalEnergy + m.VfromL(xList[x])
			// 	}
			// 	eMinCrossingIndex, eMaxCrossingIndex := int(math.Ceil(energyMin/m.Parameters.EnergyDiscretizationStep)), int(math.Ceil(energyMax/m.Parameters.EnergyDiscretizationStep))
			// 	eList := make([]float64, eMaxCrossingIndex-eMinCrossingIndex+1)
			// 	{
			// 		for e := eMinCrossingIndex; e <= eMaxCrossingIndex; e++ {
			// 			eList[e-eMinCrossingIndex] = float64(e) * m.Parameters.EnergyDiscretizationStep
			// 		}
			// 	}
			// 	xListFromE := make([]float64, len(eList))
			// 	for e := range eList {
			// 		xListFromE[e] = m.LfromV(-(tUpdate.totalEnergy - eList[e]))
			// 	}

			// 	mergedEList := make([]float64, len(xList)+len(eList)+2)
			// 	mergedXList := make([]float64, len(xList)+len(eList)+2)
			// 	{
			// 		mergedEList[0] = energyMin
			// 		mergedXList[0] = xMin
			// 		eIndex := 0
			// 		xIndex := 0
			// 		i := 0
			// 		for ; xIndex < len(xList) && eIndex < len(eList); i++ {
			// 			if eList[eIndex] < eListFromX[xIndex] {
			// 				mergedEList[i+1] = eList[eIndex]
			// 				mergedXList[i+1] = xListFromE[eIndex]
			// 				eIndex++
			// 			} else {
			// 				mergedEList[i+1] = eListFromX[xIndex]
			// 				mergedXList[i+1] = xList[xIndex]
			// 				xIndex++
			// 			}
			// 		}
			// 		for ; eIndex < len(eList); i++ {
			// 			mergedEList[i+1] = eList[eIndex]
			// 			mergedXList[i+1] = xListFromE[eIndex]
			// 			eIndex++
			// 		}
			// 		for ; xIndex < len(xList); i++ {
			// 			mergedEList[i+1] = eListFromX[xIndex]
			// 			mergedXList[i+1] = xList[xIndex]
			// 			xIndex++
			// 		}
			// 		mergedEList[len(mergedEList)-1] = energyMax
			// 		mergedXList[len(mergedXList)-1] = xMax
			// 	}
			// 	for i := range len(mergedEList) - 1 {
			// 		meanX := 0.5 * (mergedXList[i] + mergedXList[i+1])
			// 		meanE := 0.5 * (mergedEList[i] + mergedEList[i+1])
			// 		x := int(math.Floor(meanX / m.XStep))
			// 		e := int(math.Floor(meanE / m.Parameters.EnergyDiscretizationStep))
			// 		cellE := (float64(e) + 0.5) * m.Parameters.EnergyDiscretizationStep
			// 		origin := tUpdate.origin
			// 		mu := math.Sqrt((meanE - tUpdate.radialEnergy) / meanE)
			// 		time := (mergedXList[i+1] - mergedXList[i]) / utils.EV2electronVelocity(meanE-tUpdate.radialEnergy)

			// 		if 0 < x && x < m.NumCells && 0 < e && e < m.NumCellsE {
			// 			for process := range m.RateLookups {
			// 				if _, exist := DetailedCollisionAtCell[process]; exist {
			// 					update := m.RateLookups[process][e] * time / m.XStep
			// 					// update :=
			// 					DetailedCollisionAtCell[process][x][origin-m.TotalElectronsEmittedOnCathode].Add(update)
			// 				}
			// 			}
			// 			ionizations := m.ionizationCS.CrossSectionAt(cellE) * mu * utils.EV2electronVelocity(cellE) * m.Parameters.GasDensity * time / m.XStep //m.RateLookups[string(lxgata.IONIZATION)][e] * time / m.XStep
			// 			CollisionAtCell_perElectron[lxgata.IONIZATION][x][origin-m.TotalElectronsEmittedOnCathode].Add(ionizations)
			// 		}
			// 	}
			// }
			stateWg.Done()
		}()

		// var wallLossFlow chan WallLossEvent
		// if m.Parameters.Volumetric {
		// 	wallLossFlow = make(chan WallLossEvent, 10000)
		// 	stateWg.Add(1)
		// 	go func() {
		// 		for wallLossUpdate := range wallLossFlow {
		// 			if wallLossUpdate.x < len(WallLossAtCell[wallLossUpdate.x]) {
		// 				WallLossAtCell[wallLossUpdate.x][wallLossUpdate.origin] += 1. // wallLossUpdate.absVx
		// 			}

		// 		}
		// 		stateWg.Done()
		// 	}()
		// }

		ooeFlow := make(chan int, m.Parameters.NElectrons*100)
		stateWg.Add(1)
		go func() {
			for ooeEvent := range ooeFlow {
				if ooeEvent < m.NumCells {
					m.OutOfEnergyAtCell[ooeEvent]++
				}
			}
			stateWg.Done()
		}()

		for _, process := range m.Parameters.CrossSectionsData().GetTypes() {
			for x := range CollisionAtCellPerAvalanche[process] {
				CollisionAtCellPerAvalanche[process][x] = make([]int, nElectrons)
			}
		}

		// for p := range m.Parameters.CrossSectionsData().Processes {
		// 	for substring := range m.Parameters.RequireCollisionRelativeMargin {
		// 		collisionName := string(m.Parameters.CrossSectionsData().Processes[p].Type) + m.Parameters.CrossSectionsData().Processes[p].Outcome
		// 		if strings.Contains(collisionName, substring) {
		// 			for x := range DetailedCollisionAtCell[collisionName] {
		// 				DetailedCollisionAtCell[collisionName][x] = make([]utils.KahanSummable, nElectrons)
		// 			}
		// 		}
		// 	}
		// }

		// WallLossAtCell = make([]utils.WeightedDistribution, nElectrons)
		// for e := range nElectrons {
		// 	WallLossAtCell[e] = utils.NewWeightedDistribution(m.NumCells + 1)
		// }

		// for x := range TimeAtCellPerAvalanche {
		// 	for e := range TimeAtCellPerAvalanche[x] {
		// 		TimeAtCellPerAvalanche[x][e] = make([]utils.KahanSummable, nElectrons)
		// 	}
		// }

		for origin := range nElectrons {
			particle := m.newParticle(origin + m.TotalElectronsEmittedOnCathode)
			if m.Parameters.CalculateDistribution {
				// angleCell := m.getAngleCell(particle.mu)
				m.DistributionXEMu[0][int(particle.eKinetic/m.Parameters.EnergyDiscretizationStep)][int((particle.mu+1)/m.Parameters.MuDiscretizationStep)]++
				// m.Distribution[0][angleCell] = append(m.Distribution[0][int(particle.mu/m.Parameters.MuDiscretizationStep)], int(particle.eKinetic/m.EStep))
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
					// flow := make([]CellUpdate, 0)
					for lowerEnergyThreshold < particlePtr.trajectory.totEnergy { //xCell :=int((particlePtr.x+0.5*m.XStep)/m.XStep); (!m.Parameters.ParallelPlaneHollowCathode && xCell < m.NumCells) ||
						// (m.Parameters.ParallelPlaneHollowCathode && int(particlePtr.x/m.XStep) < m.NumCells+1) {
						collision, throwOut := m.nextCollision(particlePtr, trajFlow /*, stateflow*/)
						// if m.Parameters.CalculateDistribution {
						// 	flow = append(flow, distUpdate...)
						// }
						if throwOut {
							if particlePtr.x == m.Parameters.SimulationLength() {
								energyDepositionFlow <- particlePtr.trajectory.totEnergy
							}
							// if m.Parameters.Volumetric {
							// 	wallLossFlow <- WallLossEvent{int(particlePtr.x / m.XStep), utils.EV2electronVelocity(particlePtr.getAxialEnergy()), particlePtr.origin}
							// }

							break
						}
						if collision != nil {

							energyLoss := collision.Threshold
							cosChiScattered := m.Parameters.CrossSectionsData().SampleScatteringAngleCos(particlePtr.eKinetic, energyLoss, collision.Type, lxgata.IgnoreAtomicNumber)
							phi := 2. * math.Pi * rand.Float64()

							axe := particlePtr.getAxialEnergy()
							absoluteXVelocityBeforeCollision := utils.EV2electronVelocity(axe)
							if math.IsNaN(absoluteXVelocityBeforeCollision) {
								fmt.Println("avx is nan")
							}
							particlePtr.eKinetic -= energyLoss
							switch collision.Type {
							case lxgata.ELASTIC:
								energyLoss = particlePtr.eKinetic * 2. * collision.MassRatio * (1. - cosChiScattered)
								particlePtr.eKinetic -= energyLoss

							// case lxgata.EFFECTIVE:
							// 	energyLoss = particlePtr.eKinetic * 2. * collision.MassRatio * (1. - cosChiScattered)
							// 	particlePtr.eKinetic -= energyLoss

							// case lxgata.EXCITATION:
							// case lxgata.ATTACHMENT:
							case lxgata.IONIZATION:
								ejected := *particlePtr
								ejected.ejectedFromIonization = true

								availableEnergy := particlePtr.eKinetic

								switch m.Parameters.GetIonizationEnergySharingMode() {
								case config.Equal:
									ejected.eKinetic = 0.5 * particlePtr.eKinetic
								case config.Opal:
									w := utils.OpalIonizationShapeParameters[m.Parameters.Species]
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
									cosChiScattered = lxgata.CoulombScatteringAngleSample(particlePtr.eKinetic+energyLoss, m.Parameters.CrossSectionsData().UParameter, energyLoss, lxgata.IgnoreAtomicNumber)
									cosChiEjected = 1. - 2.*rand.Float64()
								case config.CoulombM:
									cosChiScattered = lxgata.CoulombScatteringAngleSample(particlePtr.eKinetic+energyLoss, m.Parameters.CrossSectionsData().UParameter, energyLoss, lxgata.IgnoreAtomicNumber)
									cosChiEjected = lxgata.CoulombScatteringAngleSample(ejected.eKinetic+energyLoss, m.Parameters.CrossSectionsData().UParameter, energyLoss, lxgata.IgnoreAtomicNumber)
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
							case lxgata.EXCITATION:
							case lxgata.ROTATION:
							default:
								panic(fmt.Sprintf("unexpected lxgata.CollisionType: %#v", collision.Type))
							}
							if !m.Parameters.Volumetric || (particlePtr.y*particlePtr.y+particlePtr.z*particlePtr.z < m.cathodeRadius2) || particlePtr.x < m.Parameters.CathodeFallLength {

								collFlow <- CollisionEvent{
									x:          int(particlePtr.x / m.XStep),
									energyLoss: energyLoss,
									collType:   collision.Type,
									origin:     particlePtr.origin,
									outcome:    collision.Outcome,
								}
							} // else {
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
						} else {
							// currentCell := int(particlePtr.x / m.XStep)
							// if currentCell < m.NumCells && m.Parameters.CountNulls {
							// 	select {
							// 	case collFlow <- CollisionEvent{currentCell, utils.EV2electronVelocity(particlePtr.getAxialEnergy()), 0, "NULL", "", particlePtr.origin}:
							// 	default:
							// 	}

							// }
						}
					}
					if particlePtr.trajectory.totEnergy < lowerEnergyThreshold {
						ooeFlow <- int(particlePtr.x / m.XStep)
					}
					// if m.Parameters.CalculateDistribution {
					// 	distFlow <- flow
					// }
					computeWg.Done()
				}
			}()
		}
		computeWg.Wait()

		close(computeFlow)
		close(collFlow)
		// close(outcollFlow)
		close(energyDepositionFlow)
		// if m.Parameters.Volumetric {
		// 	close(wallLossFlow)
		// }
		// if m.Parameters.CalculateDistribution {
		// 	close(distFlow)
		// }
		close(trajFlow)
		close(ooeFlow)
		stateWg.Wait()

		// for x := range m.TimeAtCell {
		// 	for energy := range m.TimeAtCell[x] {
		// 		m.TimeAtCell[x][energy].Update()
		// 	}
		// }

		// if m.Parameters.CellTimeWeighting {
		// 	for process := range m.CollisionAtCell {
		// 		for x := range m.NumCells {
		// 			// m.CollisionAtCell[process][x].UpdateIK(CollisionAtCellPerAvalanche[process][x], TimeAtCellPerAvalanche[x])
		// 		}
		// 		// m.CollisionOutside.GetPointer(process).Update(CollisionOutside_perElectron[process])
		// 	}
		// 	for process := range m.DetailedCollisionAtCell {
		// 		for x := range m.NumCells {
		// 			// m.DetailedCollisionAtCell[process][x].UpdateIK(DetailedCollisionAtCell[process][x], TimeAtCellPerAvalanche[x])
		// 		}
		// 	}

		// 	for x := range m.TimeAtCell {
		// 		m.TimeAtCell[x].UpdateK(TimeAtCellPerAvalanche[x])
		// 	}
		// } else {
		for process := range m.CollisionAtCell {
			for x := range m.NumCells {
				m.CollisionAtCell[process][x].Update(CollisionAtCellPerAvalanche[process][x])
			}
			// m.CollisionOutside.GetPointer(process).Update(CollisionOutside_perElectron[process])
		}
		for process := range m.DetailedCollisionAtCell {
			for x := range m.NumCells {
				m.DetailedCollisionAtCell[process][x].Update(DetailedCollisionAtCell[process][x])
			}
		}
		// }

		// m.WallLossAtCell.Update(WallLossAtCell)

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
			// if m.Parameters.CellTimeWeighting {
			// 	m.IonizationsSumUpToCell[x].UpdateIK(PerAvalancheCollisionSumsUpToCell, PerAvalancheTimeSumsUpToCell)
			// } else {
			m.IonizationsSumUpToCell[x].Update(PerAvalancheCollisionSumsUpToCell)
			// }

		}
		m.TotalElectronsEmittedOnCathode += nElectrons
		nElectrons = electronsToSimulate(m)
	}

	if numberOfEDEvents > 0 {
		m.MeanElectronEnergyAtAnode = energyDeposition / float64(numberOfEDEvents)
	}
	print("\r")
}
