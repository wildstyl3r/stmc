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
	// boundaryThickness        float64
	cathodeRadius2, tubeRadius2 float64
	// boundaryRadius2          float64

	Va float64 // additional voltage to avoid numeric negative energy beyond cathode fall

	NumCells, NumCellsMu, NumCellsE int

	XStep float64

	// lookUpPotential         []float64
	LookupEnergy []float64
	// LookupInverseVelocity     []float64
	// GridInverseMu           []float64
	lookupMNumerator []float64

	Distribution [][][]int

	CollisionAtCell         map[lxgata.CollisionType][]utils.Aggregation
	CollisionOutside        map[lxgata.CollisionType][]utils.Aggregation
	DetailedCollisionAtCell map[string][]utils.Aggregation
	MeanEnergy              []utils.Aggregation
	IonizationsSumUpToCell  []utils.Aggregation
	WallLossAtCell          []utils.Aggregation
	EnergyLossByProcess     map[lxgata.CollisionType][]float64

	lookupCosinesAtAngleCellBounds  []float64
	LookupCosinesAtAngleCellCenters []float64

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
	m.constFieldAcceleration = constants.ElementaryCharge * parameters.ConstEField / constants.ElectornMass

	m.sheathTimeOuterConstant = 2 * parameters.CathodeFallLength * math.Sqrt(constants.ElectornMass/(2*constants.ElementaryCharge*m.Vc))
	m.sheathTimeInnerConstant = math.Sqrt(parameters.CathodeFallPotential) / parameters.CathodeFallLength
	m.constFieldTimeConstant = math.Sqrt(2*constants.ElectornMass/constants.ElementaryCharge) * (-parameters.ConstEField)

	m.CalculateStaticTimeMotionConstants()

	meanFreePath := 1. / (m.Parameters.CrossSectionsData().SurplusCrossSection() * m.Parameters.GasDensity)
	// if m.Parameters.Verbose() {
	// 	fmt.Printf("Mean free path: %f\n", meanFreePath)
	// }

	m.NumCells = 5 * int(m.Parameters.SimulationLength()/meanFreePath)
	m.NumCellsMu = int(2./parameters.MuDiscretizationStep) + 1
	m.NumCellsE = int((m.Va + m.Vc + 5) / m.Parameters.EnergyDiscretizationStep)
	// if m.Parameters.Verbose() {
	// 	fmt.Printf("xCells: %d;\n", m.NumCells)
	// }

	m.XStep = m.Parameters.SimulationLength() / float64(m.NumCells)
	// m.MuStep = m.Parameters.MuStep    // eV
	// m.NumMuCells = int(2. / m.MuStep)

	radianAngleStep := parameters.AngleStep / 360. * 2. * math.Pi
	numAngleCells := int(360. / parameters.AngleStep)
	m.lookupCosinesAtAngleCellBounds, m.LookupCosinesAtAngleCellCenters = make([]float64, numAngleCells), make([]float64, numAngleCells)
	for a := range numAngleCells {
		angle := float64(a) * radianAngleStep
		nextAngle := float64(a+1) * radianAngleStep
		m.lookupCosinesAtAngleCellBounds[a] = math.Cos(angle)
		m.LookupCosinesAtAngleCellCenters[a] = math.Cos(0.5 * (angle + nextAngle))
	}

	m.CollisionAtCell = make(map[lxgata.CollisionType][]utils.Aggregation)
	m.CollisionOutside = make(map[lxgata.CollisionType][]utils.Aggregation)
	m.DetailedCollisionAtCell = make(map[string][]utils.Aggregation)
	m.EnergyLossByProcess = make(map[lxgata.CollisionType][]float64)
	m.IonizationsSumUpToCell = make([]utils.Aggregation, m.NumCells+1)
	m.MeanEnergy = make([]utils.Aggregation, m.NumCells+1)
	for _, process := range parameters.CrossSectionsData().GetTypes() {
		m.EnergyLossByProcess[process] = make([]float64, m.NumCells+1)
		m.CollisionAtCell[process] = make([]utils.Aggregation, m.NumCells+1)
		m.CollisionOutside[process] = make([]utils.Aggregation, m.NumCells+1)
	}
	m.DetailedCollisionAtCell = make(map[string][]utils.Aggregation)
	for p := range parameters.CrossSectionsData().Processes {
		for substring := range parameters.RequireCollisionRelativeMargin {
			collisionName := string(parameters.CrossSectionsData().Processes[p].Type) + parameters.CrossSectionsData().Processes[p].Outcome
			if strings.Contains(collisionName, substring) {
				m.DetailedCollisionAtCell[collisionName] = make([]utils.Aggregation, m.NumCells+1)
			}
		}

	}
	m.WallLossAtCell = make([]utils.Aggregation, m.NumCells+1)

	// m.lookUpPotential = make([]float64, m.NumCells+2)
	// for cellNode := range m.lookUpPotential {
	// m.lookUpPotential[cellNode] = m.VfromL(float64(cellNode) * m.XStep, false)
	// }

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
		m.Distribution = make([][][]int, m.NumCells)
		for x := range m.Distribution {
			m.Distribution[x] = make([][]int, m.NumCellsE)
			for e := range m.Distribution[x] {
				m.Distribution[x][e] = make([]int, m.NumCellsMu)
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

func (m *Model) SheathTime(x1, axialEnergy1_eV, axialEnergy2_eV float64) float64 { //time to get from x1 to x2 assuming both are in the sheath
	c1, c2 := m.GetDynamicTimeMotionConstants(x1, axialEnergy1_eV, (axialEnergy2_eV - axialEnergy1_eV))
	maxWT := utils.TernarySearchMax(func(wt float64) float64 {
		sht := m.SheathEnergyAtWTime(wt, c1, c2)
		delta := sht - axialEnergy2_eV
		return delta
	}, 0, math.Pi, 1e-9)
	Wt, _ := utils.Bisect(func(wt float64) float64 {
		sht := m.SheathEnergyAtWTime(wt, c1, c2)
		delta := sht - axialEnergy2_eV
		return delta
	}, 0, maxWT, 1e-9, 1e-9)
	t := Wt / m.timeMotionW
	if -t > 1e-9 {
		fmt.Printf("sheath t < 0\n")
	}
	return t
}

func (m *Model) ConstEFieldTime(axialEnergy1_eV, axialEnergy2_eV float64) float64 { //time to get from x1 to x2 in the constant field region
	return m.constFieldTimeConstant * math.Abs(math.Sqrt(axialEnergy2_eV)-math.Sqrt(axialEnergy1_eV))
}

// func (m *Model) getAngleCell(mu float64) int {
// 	for a := range m.lookupCosinesAtAngleCellBounds {
// 		if a+1 < len(m.lookupCosinesAtAngleCellBounds) && m.lookupCosinesAtAngleCellBounds[a] < mu && mu < m.lookupCosinesAtAngleCellBounds[a+1] {
// 			return a
// 		}
// 	}
// 	return len(m.lookupCosinesAtAngleCellBounds) - 1
// }

func (m *Model) collisionSelector(eKinetic, x, M float64) *lxgata.Collision {
	// var crossSections = m.Parameters.CrossSectionsData().CrossSectionsAt(eKinetic)
	var totalCrossSectionPrimed = M * math.Abs(m.EFieldFromL(x)) / math.Sqrt(eKinetic) / m.Parameters.GasDensity
	return m.Parameters.CrossSectionsData().SampleWithNullCollision(eKinetic, totalCrossSectionPrimed)
	// var crossSectionAccum float64 = 0. //totalCrossSectionPrimed - totalCrossSection // the difference is null collision
	// var choice float64 = rand.Float64() * totalCrossSectionPrimed

	// for i := range m.Parameters.CrossSectionsData().Processes {
	// 	crossSectionAccum += crossSections[i]
	// 	if choice < crossSectionAccum {
	// 		return &(m.Parameters.CrossSectionsData().Processes[i])
	// 	}
	// }
	// return nil
}

type FlowElem struct {
	x, mu, energy int
}

type Flow []FlowElem

func (m *Model) getFlowBetweenEnergies(startEnergy, endEnergy, totalEnergy, radialEnergy float64) (flow []FlowElem) {
	energyMin, energyMax := min(startEnergy, endEnergy), max(startEnergy, endEnergy)
	var direction float64
	if startEnergy < endEnergy {
		direction = +1.
	} else {
		direction = -1.
	}
	xMin, xMax := m.LfromV(-(totalEnergy - energyMin)), m.LfromV(-(totalEnergy - energyMax))
	gridMinIndex, gridMaxIndex := int(xMin/m.XStep)+1, int(xMax/m.XStep)+1
	if gridMinIndex < 0 {
		panic("i<0")
	}
	for i := gridMinIndex; i < gridMaxIndex; i++ {
		eKinetic := totalEnergy + m.VfromL(float64(i)*m.XStep)
		mu := math.Copysign(math.Sqrt((eKinetic-radialEnergy)/eKinetic), direction)
		flow = append(flow, FlowElem{i, int((mu + 1) / m.Parameters.MuDiscretizationStep), int(eKinetic / m.Parameters.EnergyDiscretizationStep)})
	}
	return
}

func (m *Model) timeToWall(p *Particle) float64 {
	if math.IsInf(m.Parameters.TubeRadius, 0) {
		return math.Inf(-1)
	}
	vR := utils.EV2electronVelocity(p.eStar)
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

func (m *Model) GetDynamicTimeMotionConstants(x0, ex0, c2sign float64) (C1, C2 float64) {
	C1 = x0 - m.timeMotionP
	C2 = math.Copysign(m.Parameters.CathodeFallLength*math.Sqrt(ex0/m.timeMotionK), c2sign)
	return
}

func (m *Model) SheathPositionAtTime(t, C1, C2 float64) (x float64) {
	return C1*math.Cos(m.timeMotionW*t) + C2*math.Sin(m.timeMotionW*t) + m.timeMotionP
}

func (m *Model) SheathEnergyAtWTime(wt, C1, C2 float64) (e float64) {
	sp := C2*math.Cos(wt) - C1*math.Sin(wt)
	energyEV := m.timeMotionK * sp * sp * m.inverseCathodeFallLength * m.inverseCathodeFallLength
	return energyEV
}

func (m *Model) WallIntersection(p *Particle) (tw, xw float64) {
	C1, C2 := m.GetDynamicTimeMotionConstants(p.x, p.getAxialEnergy(), p.mu)
	tw = m.timeToWall(p) * m.timeMotionW
	if math.IsInf(tw, 0) {
		return
	}
	xw = m.SheathPositionAtTime(tw, C1, C2)
	if xw < m.Parameters.CathodeFallLength {
		return
	}
	axEnergyAtDc := p.getAxialEnergy() - (m.VfromL(m.Parameters.CathodeFallLength) - m.VfromL(p.x))
	td := m.SheathTime(p.x, p.getAxialEnergy(), axEnergyAtDc)
	tc := tw - td
	xw = m.Parameters.CathodeFallLength + utils.EV2electronVelocity(axEnergyAtDc)*tc + 0.5*m.constFieldAcceleration*tc*tc
	return
}

func (m *Model) nextCollision(p *Particle) (collisionDescription *lxgata.Collision, flow []FlowElem, throwOut bool) {
	wallCollisionEnergyStep := -1
	var tw, xw, ew float64
	if m.Parameters.Volumetric {
		tw, xw = m.WallIntersection(p)
		if !math.IsInf(tw, 0) {
			ew = p.totEnergy + m.VfromL(xw)
			wallCollisionEnergyStep = int(ew / m.Parameters.EnergyStep)
		}
	}

	R := -math.Log(1. - rand.Float64())

	if !(0 <= p.x && p.x < m.Parameters.SimulationLength()) {
		throwOut = true
		return
	}

	var currentCellIndex int
	var cachedVelocity float64
	var isVelocityCached = false

	alignedToEnergyGrid := false

	for {
		var lowEnergy, segmentStartEnergy, segmentEndEnergy, highEnergy float64
		var nextCellIndex, highEnergyCellIndex int
		var reversalOccured, arrivalAtCathode, arrivalAtGapEnd, arrivalAtTubeWall, highEnergyAligned bool
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
			if lowEnergy < p.eStar {
				lowEnergy, alignedToEnergyGrid = p.eStar, false
				nextCellIndex = currentCellIndex
				reversalOccured = true
			}

			//check for cathode return
			if p.totEnergy-lowEnergy > m.Vc+m.Va {
				arrivalAtCathode = true
				lowEnergy = p.totEnergy - (m.Vc + m.Va)
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
			if p.totEnergy <= highEnergy {
				arrivalAtGapEnd = true
				highEnergy, highEnergyAligned = p.totEnergy, false
			} else {
				highEnergyCellIndex, highEnergyAligned = nextCellIndex, true
			}
			segmentStartEnergy, segmentEndEnergy = lowEnergy, highEnergy
		}

		var M, G float64
		if highEnergyAligned {
			M = m.lookupMNumerator[highEnergyCellIndex] / -m.EFieldFromPotential(-(p.totEnergy - highEnergy))
		} else {
			M = m.Parameters.GasDensity * m.Parameters.CrossSectionsData().TotalCrossSectionAt(highEnergy) * math.Sqrt(highEnergy) / -m.EFieldFromPotential(-(p.totEnergy - highEnergy))
		}

		var higherVelocity, lowerVelocity float64
		if isVelocityCached {
			if p.mu < 0 {
				higherVelocity, lowerVelocity = cachedVelocity, math.Sqrt(lowEnergy-p.eStar)
			} else {
				higherVelocity, lowerVelocity = math.Sqrt(highEnergy-p.eStar), cachedVelocity
			}
		} else {
			higherVelocity, lowerVelocity = math.Sqrt(highEnergy-p.eStar), math.Sqrt(lowEnergy-p.eStar)
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
		} else {
			{
				var delta float64
				if p.mu < 0 {
					delta = math.Sqrt(segmentStartEnergy-p.eStar) - R/(2.*M) // == sqrt(pColl - p.eStar) i.e. abs(speed) equivalent along x axis
					if delta < 0 {
						R -= 2 * M * math.Sqrt(segmentStartEnergy-p.eStar)
						segmentEndEnergy = p.eStar
						reversalOccured = true
					}
					if math.IsNaN(delta) {
						panic("delta is NaN")
					}
				} else {
					delta = R/(2.*M) + math.Sqrt(segmentStartEnergy-p.eStar)
					if math.IsNaN(delta) {
						panic("delta is NaN")
					}
				}
				if !reversalOccured {
					collisionEnergy := math.FMA(delta, delta, p.eStar) // R = 2M[sqrt(p.e-p.eStar) - sqrt(eColl - p.eStar)]
					segmentEndEnergy = collisionEnergy
					p.setEnergy(collisionEnergy, m, true, true)
					if p.totEnergy < collisionEnergy || m.Parameters.SimulationLength() < p.x {
						arrivalAtGapEnd = true
						segmentEndEnergy = p.totEnergy
					} else if p.x < 0 {
						arrivalAtCathode = true
						segmentEndEnergy = p.totEnergy - (m.Vc + m.Va)
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

		if m.Parameters.CalculateDistribution {
			flow = append(flow, m.getFlowBetweenEnergies(segmentStartEnergy, segmentEndEnergy, p.totEnergy, p.eStar)...)
		}

		if collisionOccured {
			return
		}

		if arrivalAtCathode {
			p.x = 0
			throwOut = true
			return
		}

		if arrivalAtTubeWall {
			p.x = xw
			var tubeBackscatter0 float64 //tubeBackscatterB
			if m.Parameters.ParallelPlaneHollowCathode {
				tubeBackscatter0 = m.Parameters.AnodeBackscatteringCoefficient0
				// tubeBackscatterB = m.Parameters.AnodeBackscatteringCoefficientB
			} else {
				tubeBackscatter0 = m.Parameters.TubeBackscatteringCoefficient0
				// tubeBackscatterB = m.Parameters.TubeBackscatteringCoefficientB
			}
			if tubeBackscatter0 != 0 {
				// backscatteringCoefficient := tubeBackscatter0 * math.Exp(tubeBackscatterB*(1-math.Cos(incid)))
			} else {
				throwOut = true
				return
			}
		}

		if arrivalAtGapEnd {
			if m.Parameters.ParallelPlaneHollowCathode {
				if p.mu < 0 {
					println("DEBUG: arrival at half-gap from beyond")
				}
				p.mu = -p.mu
				// p.eKinetic = p.totEnergy
				p.setEnergy(p.totEnergy, m, true, false)
				alignedToEnergyGrid = false
			} else {
				if m.Parameters.AnodeBackscatteringCoefficient0 != 0 {
					backscatteringCoefficient := m.Parameters.AnodeBackscatteringCoefficient0 * math.Exp(m.Parameters.AnodeBackscatteringCoefficientB*(1-p.mu))

					p.setEnergy(p.totEnergy, m, true, true)
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
					p.eKinetic = p.totEnergy
					p.mu = math.Sqrt(p.getAxialEnergy() / p.totEnergy)
					return
				}

			}
		} else if reversalOccured {
			p.mu = +0.
			p.eKinetic = p.eStar
			alignedToEnergyGrid = false
		} else {
			alignedToEnergyGrid = true
		}
		if m.Parameters.Volumetric && currentCellIndex == wallCollisionEnergyStep && !math.IsInf(m.tubeRadius2, 0) {
			throwOut = true
			return
		}
		currentCellIndex = nextCellIndex
	}
}

type CollisionEvent struct {
	x          int
	absVx      float64
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
	CollisionAtCell := make(map[lxgata.CollisionType][][]int)
	CollisionOutside := make(map[lxgata.CollisionType][][]int)
	for _, process := range m.Parameters.CrossSectionsData().GetTypes() {
		CollisionAtCell[process] = make([][]int, m.NumCells+1)
		CollisionOutside[process] = make([][]int, m.NumCells+1)
	}
	DetailedCollisionAtCell := make(map[string][][]int)
	for p := range m.Parameters.CrossSectionsData().Processes {
		for substring := range m.Parameters.RequireCollisionRelativeMargin {
			collisionName := string(m.Parameters.CrossSectionsData().Processes[p].Type) + m.Parameters.CrossSectionsData().Processes[p].Outcome
			if strings.Contains(collisionName, substring) {
				DetailedCollisionAtCell[collisionName] = make([][]int, m.NumCells+1)
			}
		}
	}
	WallLossAtCell := make([][]int, m.NumCells+1)
	IonizationsAtCellE := make([][]int, m.NumCells+1)

	var energyDeposition float64
	var numberOfEDEvents int

	nElectrons := electronsToSimulate(m)
	for nElectrons > 0 {
		m.PurgeMetrics()

		computeFlow := make(chan *Particle, nElectrons*100)
		var computeWg, stateWg sync.WaitGroup

		collFlow := make(chan CollisionEvent, 100000)
		stateWg.Add(1)
		go func() {
			for collision := range collFlow {
				if collision.x < m.NumCells+1 {
					CollisionAtCell[collision.collType][collision.x][collision.origin-m.TotalElectronsEmittedOnCathode]++ // collision.absVx

					if collision.collType == lxgata.IONIZATION {
						IonizationsAtCellE[collision.x][collision.origin-m.TotalElectronsEmittedOnCathode]++ // collision.absVx
					}

					if collision.collType != "NULL" {
						m.EnergyLossByProcess[collision.collType][collision.x] += collision.energyLoss
						detailedName := string(collision.collType) + collision.outcome
						if _, exist := m.DetailedCollisionAtCell[detailedName]; exist {
							DetailedCollisionAtCell[detailedName][collision.x][collision.origin-m.TotalElectronsEmittedOnCathode]++ // collision.absVx
						}
					}
				}
			}
			stateWg.Done()
		}()

		outcollFlow := make(chan CollisionEvent, 100000)
		stateWg.Add(1)
		go func() {
			for collision := range outcollFlow {
				if collision.x < m.NumCells+1 {
					CollisionOutside[collision.collType][collision.x][collision.origin-m.TotalElectronsEmittedOnCathode]++ // collision.absVx
				}
			}
			stateWg.Done()
		}()

		energyDepositionFlow := make(chan float64, 10*m.Parameters.NElectrons)
		stateWg.Add(1)
		go func() {
			for deposition := range energyDepositionFlow {
				energyDeposition += deposition
				numberOfEDEvents += 1
			}
			stateWg.Done()
		}()

		var distFlow chan Flow
		if m.Parameters.CalculateDistribution {
			distFlow = make(chan Flow, 10000)
			stateWg.Add(1)
			go func() {
				for distUpdate := range distFlow {
					for i := range distUpdate {
						x, e, mu := distUpdate[i].x, distUpdate[i].energy, distUpdate[i].mu
						if x < m.NumCells && mu < m.NumCellsMu && e < m.NumCellsE {
							m.Distribution[x][e][mu]++
						}
					}
				}
				stateWg.Done()
			}()
		}

		var wallLossFlow chan WallLossEvent
		if m.Parameters.Volumetric {
			wallLossFlow = make(chan WallLossEvent, 10000)
			stateWg.Add(1)
			go func() {
				for wallLossUpdate := range wallLossFlow {
					if wallLossUpdate.x < len(WallLossAtCell[wallLossUpdate.x]) {
						WallLossAtCell[wallLossUpdate.x][wallLossUpdate.origin] += 1. // wallLossUpdate.absVx
					}

				}
				stateWg.Done()
			}()
		}

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
			for c := range m.NumCells + 1 {
				CollisionAtCell[process][c] = make([]int, nElectrons)
				CollisionOutside[process][c] = make([]int, nElectrons)
			}
		}
		for collisionName := range DetailedCollisionAtCell {
			for c := range m.NumCells + 1 {
				DetailedCollisionAtCell[collisionName][c] = make([]int, nElectrons)
			}
		}

		for c := range m.WallLossAtCell {
			WallLossAtCell[c] = make([]int, nElectrons)
			IonizationsAtCellE[c] = make([]int, nElectrons)
		}

		for origin := range nElectrons {
			particle := m.newParticle(origin + m.TotalElectronsEmittedOnCathode)
			if m.Parameters.CalculateDistribution {
				// angleCell := m.getAngleCell(particle.mu)
				m.Distribution[0][int(particle.eKinetic/m.Parameters.EnergyDiscretizationStep)][int((particle.mu+1)/m.Parameters.MuDiscretizationStep)]++
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
					flow := make([]FlowElem, 0)
					for lowerEnergyThreshold < particlePtr.totEnergy || m.Parameters.CalculateDistribution { //xCell :=int((particlePtr.x+0.5*m.XStep)/m.XStep); (!m.Parameters.ParallelPlaneHollowCathode && xCell < m.NumCells) ||
						// (m.Parameters.ParallelPlaneHollowCathode && int(particlePtr.x/m.XStep) < m.NumCells+1) {
						collision, distUpdate, throwOut := m.nextCollision(particlePtr /*, stateflow*/)
						if m.Parameters.CalculateDistribution {
							flow = append(flow, distUpdate...)
						}
						if throwOut {
							if particlePtr.x == m.Parameters.SimulationLength() {
								energyDepositionFlow <- particlePtr.totEnergy
							}
							if m.Parameters.Volumetric {
								wallLossFlow <- WallLossEvent{int(particlePtr.x / m.XStep), utils.EV2electronVelocity(particlePtr.getAxialEnergy()), particlePtr.origin}
							}

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
								ejected._debug_IonEjected = true

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
								if ejected.totEnergy >= lowerEnergyThreshold || m.Parameters.CalculateDistribution {
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
									absVx:      absoluteXVelocityBeforeCollision,
									energyLoss: energyLoss,
									collType:   collision.Type,
									origin:     particlePtr.origin,
									outcome:    collision.Outcome,
								}
							} else {
								outcollFlow <- CollisionEvent{
									x:          int(particlePtr.x / m.XStep),
									absVx:      absoluteXVelocityBeforeCollision,
									energyLoss: energyLoss,
									collType:   collision.Type,
									origin:     particlePtr.origin,
									outcome:    collision.Outcome,
								}
							}

							if collision.Type == lxgata.ATTACHMENT {
								break
							}
							particlePtr.redirect(cosChiScattered, math.Cos(phi), m)
						} else {
							currentCell := int(particlePtr.x / m.XStep)
							if currentCell < m.NumCells && m.Parameters.CountNulls {
								select {
								case collFlow <- CollisionEvent{currentCell, utils.EV2electronVelocity(particlePtr.getAxialEnergy()), 0, "NULL", "", particlePtr.origin}:
								default:
								}

							}
						}
					}
					if particlePtr.totEnergy < lowerEnergyThreshold {
						ooeFlow <- int(particlePtr.x / m.XStep)
					}
					if m.Parameters.CalculateDistribution {
						distFlow <- flow
					}
					computeWg.Done()
				}
			}()
		}
		computeWg.Wait()

		close(computeFlow)
		close(collFlow)
		close(outcollFlow)
		close(energyDepositionFlow)
		if m.Parameters.Volumetric {
			close(wallLossFlow)
		}
		if m.Parameters.CalculateDistribution {
			close(distFlow)
		}
		close(ooeFlow)
		stateWg.Wait()

		for process := range m.CollisionAtCell {
			for c := range m.NumCells + 1 {
				m.CollisionAtCell[process][c].Update(CollisionAtCell[process][c])
				m.CollisionOutside[process][c].Update(CollisionOutside[process][c])
			}
		}
		for process := range m.DetailedCollisionAtCell {
			for c := range m.NumCells + 1 {
				m.DetailedCollisionAtCell[process][c].Update(DetailedCollisionAtCell[process][c])
			}
		}
		for c := range m.WallLossAtCell {
			m.WallLossAtCell[c].Update(WallLossAtCell[c])
		}

		PerElectronSumsUpToCell := make([]int, nElectrons)
		for c := range m.IonizationsSumUpToCell {
			for e := range nElectrons {
				PerElectronSumsUpToCell[e] += IonizationsAtCellE[c][e]
			}
			m.IonizationsSumUpToCell[c].Update(PerElectronSumsUpToCell)
		}
		m.TotalElectronsEmittedOnCathode += nElectrons
		nElectrons = electronsToSimulate(m)
	}

	if numberOfEDEvents > 0 {
		m.MeanElectronEnergyAtAnode = energyDeposition / float64(numberOfEDEvents)
	}
	print("\r")
}
