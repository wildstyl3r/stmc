package model

import (
	"fmt"
	"math"
	"math/rand"
	"sync"

	"github.com/wildstyl3r/lxgata"
	"github.com/wildstyl3r/stmc/internal/config"
	"github.com/wildstyl3r/stmc/internal/utils"
)

type Model struct {
	Parameters               config.ModelParameters
	Vc                       float64 // cathode fall potential
	inverseAbsConstEField    float64
	inverseNegativeVc        float64
	inverseCathodeFallLength float64

	Va float64 // additional voltage to avoid numeric negative energy beyond cathode fall

	NumCells, NumCellsMu, NumCellsE int

	XStep float64

	lookUpPotential         []float64
	lookupTotalCrossSection []float64 // CS[energy]
	LookupEnergy            []float64
	LookupInverseVelocity   []float64
	// GridInverseMu           []float64
	lookupMNumerator []float64

	Distribution [][][]int

	CollisionAtCell     map[lxgata.CollisionType][][]uint16
	EnergyLossByProcess map[lxgata.CollisionType][]float64

	lookupCosinesAtAngleCellBounds  []float64
	LookupCosinesAtAngleCellCenters []float64

	OutOfEnergyAtCell []int
	DataHub           *DataHubType
}

func NewModel(parameters config.ModelParameters) Model {
	m := Model{}
	if parameters.ParallelPlaneHollowCathode {
		parameters.GapLength = parameters.GapLength / 2
	}
	m.Parameters = parameters
	m.Vc = math.Abs(m.Parameters.CathodeFallPotential)
	m.Va = math.Abs(m.Parameters.ConstEField * (m.Parameters.GapLength - m.Parameters.CathodeFallLength))
	m.inverseAbsConstEField = math.Abs(1. / m.Parameters.ConstEField)
	m.inverseNegativeVc = -1. / m.Vc
	m.inverseCathodeFallLength = 1. / m.Parameters.CathodeFallLength

	meanFreePath := 1. / (m.Parameters.CrossSectionsData().SurplusCrossSection() * m.Parameters.GasDensity)
	if m.Parameters.Verbose() {
		fmt.Printf("Mean free path: %f\n", meanFreePath)
	}

	m.NumCells = 5 * int(m.Parameters.GapLength/meanFreePath)
	m.NumCellsMu = int(2./parameters.MuDiscretizationStep) + 1
	m.NumCellsE = int((m.Va + m.Vc + 5) / m.Parameters.EnergyDiscretizationStep)
	if m.Parameters.Verbose() {
		fmt.Printf("xCells: %d;\n", m.NumCells)
	}

	m.XStep = m.Parameters.GapLength / float64(m.NumCells)
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

	m.CollisionAtCell = make(map[lxgata.CollisionType][][]uint16)
	m.EnergyLossByProcess = make(map[lxgata.CollisionType][]float64)
	for _, process := range parameters.CrossSectionsData().GetTypes() {
		if m.Parameters.CalculateStdError {
			m.EnergyLossByProcess[process] = make([]float64, m.NumCells+1)
			m.CollisionAtCell[process] = make([][]uint16, m.NumCells+1)
			for c := range m.NumCells + 1 {
				m.CollisionAtCell[process][c] = make([]uint16, parameters.NElectrons)
			}
		} else {
			m.EnergyLossByProcess[process] = make([]float64, m.NumCells+1)
			m.CollisionAtCell[process] = make([][]uint16, m.NumCells+1)
			for c := range m.NumCells + 1 {
				m.CollisionAtCell[process][c] = make([]uint16, 1)
			}
		}
	}
	m.lookUpPotential = make([]float64, m.NumCells+2)
	for cellNode := range m.lookUpPotential {
		m.lookUpPotential[cellNode] = m.VfromL(float64(cellNode) * m.XStep)
	}

	m.OutOfEnergyAtCell = make([]int, m.NumCells)

	numEnergyCells := int((m.Vc + m.Va + 5 + 0.1) / m.Parameters.EnergyStep)
	m.lookupTotalCrossSection = make([]float64, numEnergyCells)
	m.LookupEnergy = make([]float64, numEnergyCells)
	m.LookupInverseVelocity = make([]float64, numEnergyCells)
	m.lookupMNumerator = make([]float64, numEnergyCells)
	for gridNode := range m.LookupEnergy {
		m.LookupEnergy[gridNode] = float64(gridNode) * m.Parameters.EnergyStep
		m.lookupTotalCrossSection[gridNode] = m.Parameters.CrossSectionsData().TotalCrossSectionAt(m.LookupEnergy[gridNode])
		m.LookupInverseVelocity[gridNode] = 1. / utils.EV2electronVelocity(m.LookupEnergy[gridNode]+0.5*parameters.EnergyStep)
		m.lookupMNumerator[gridNode] = m.Parameters.GasDensity * m.lookupTotalCrossSection[gridNode] * math.Sqrt(m.LookupEnergy[gridNode])
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

	return m
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
	var crossSections = m.Parameters.CrossSectionsData().CrossSectionsAt(eKinetic)
	var totalCrossSection, totalCrossSectionPrimed = utils.SumFloat64Slice(crossSections), M * math.Abs(m.EFieldFromL(x)) / math.Sqrt(eKinetic) / m.Parameters.GasDensity
	var crossSectionAccum float64 = totalCrossSectionPrimed - totalCrossSection // the difference is null collision
	var choice float64 = rand.Float64() * totalCrossSectionPrimed

	selected := false
	if choice < crossSectionAccum {
		return nil
	}
	for i := range m.Parameters.CrossSectionsData().Processes {
		crossSectionAccum += crossSections[i]
		if choice < crossSectionAccum && !selected {
			return &(m.Parameters.CrossSectionsData().Processes[i])
		}
	}
	return nil
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

func (m *Model) nextCollision(p *Particle) (collisionDescription *lxgata.Collision, flow []FlowElem, throwOut bool) {
	R := -math.Log(1. - rand.Float64())

	if !(0 <= p.x && p.x < m.Parameters.GapLength) {
		throwOut = true
		return
	}

	var currentCellIndex int
	var cachedVelocity float64
	var isVelocityCached = false
	// potential at dc = -Va

	alignedToEnergyGrid := false

	for {
		var lowEnergy, segmentStartEnergy, segmentEndEnergy, highEnergy float64
		var nextCellIndex, highEnergyCellIndex int
		var reversalOccured, arrivalAtCathode, arrivalAtGapEnd, highEnergyAligned bool
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
				lowEnergy, highEnergy = m.LookupEnergy[currentCellIndex], m.LookupEnergy[currentCellIndex+1]
				nextCellIndex = currentCellIndex + 1
			}

			//check for gap end arrival
			if p.totEnergy < highEnergy {
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
					if p.totEnergy < collisionEnergy || m.Parameters.GapLength < p.x {
						arrivalAtGapEnd = true
						segmentEndEnergy = p.totEnergy
					}
					if p.x < 0 {
						arrivalAtCathode = true
						segmentEndEnergy = p.totEnergy - (m.Vc + m.Va)
					}
					collisionDescription = m.collisionSelector(collisionEnergy, p.x, M)
					collisionOccured = true
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
			throwOut = true
			return
		}

		if arrivalAtGapEnd {
			if m.Parameters.ParallelPlaneHollowCathode {
				if p.mu < 0 {
					println("DEBUG: arrival at half-gap from beyond")
				}
				p.mu = math.Copysign(p.mu, -1.)
				p.eKinetic = p.totEnergy
				alignedToEnergyGrid = false
			} else {
				throwOut = true
				return
			}
		}

		if reversalOccured {
			p.mu = +0.
			p.eKinetic = p.eStar
			alignedToEnergyGrid = false
		} else {
			alignedToEnergyGrid = true
		}
		currentCellIndex = nextCellIndex
	}
}

type CollisionEvent struct {
	x          int
	energyLoss float64
	collType   lxgata.CollisionType
	origin     int
}

func (m *Model) Run() {
	var computeWg, stateWg sync.WaitGroup

	collflow := make(chan CollisionEvent, 100000)
	stateWg.Add(1)
	go func() {
		for collision := range collflow {
			if collision.x < m.NumCells+1 {
				m.CollisionAtCell[collision.collType][collision.x][collision.origin]++
				if collision.collType != "NULL" {
					m.EnergyLossByProcess[collision.collType][collision.x] += collision.energyLoss
				}
			}
		}
		stateWg.Done()
	}()

	var distflow chan Flow
	if m.Parameters.CalculateDistribution {
		distflow = make(chan Flow, 10000)
		stateWg.Add(1)
		go func() {
			for distUpdate := range distflow {
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

	ooeflow := make(chan int, m.Parameters.NElectrons*100)
	stateWg.Add(1)
	go func() {
		for ooeEvent := range ooeflow {
			if ooeEvent < m.NumCells {
				m.OutOfEnergyAtCell[ooeEvent]++
			}
		}
		stateWg.Done()
	}()

	computeflow := make(chan *Particle, m.Parameters.NElectrons*10000)
	for i := range m.Parameters.NElectrons {
		origin := 0
		if m.Parameters.CalculateStdError {
			origin = i
		}
		particle := m.newParticle(origin)
		if m.Parameters.CalculateDistribution {
			// angleCell := m.getAngleCell(particle.mu)
			m.Distribution[0][int(particle.eKinetic/m.Parameters.EnergyDiscretizationStep)][int((particle.mu+1)/m.Parameters.MuDiscretizationStep)]++
			// m.Distribution[0][angleCell] = append(m.Distribution[0][int(particle.mu/m.Parameters.MuDiscretizationStep)], int(particle.eKinetic/m.EStep))
		}
		computeWg.Add(1)
		computeflow <- &particle
	}

	status := []string{"//", "==", "\\\\", "||"}
	for range m.Parameters.Threads() {
		go func() {
			counter := 0
			for particlePtr := range computeflow {
				counter++
				print("\r" + status[counter&0b11])
				ionizationThreshold := m.Parameters.CrossSectionsData().MinThresholdOfKind(lxgata.IONIZATION)
				flow := make([]FlowElem, 0)
				for ionizationThreshold < particlePtr.totEnergy || m.Parameters.CalculateDistribution { //xCell :=int((particlePtr.x+0.5*m.XStep)/m.XStep); (!m.Parameters.ParallelPlaneHollowCathode && xCell < m.NumCells) ||
					// (m.Parameters.ParallelPlaneHollowCathode && int(particlePtr.x/m.XStep) < m.NumCells+1) {
					collision, distUpdate, throwOut := m.nextCollision(particlePtr /*, stateflow*/)
					if m.Parameters.CalculateDistribution {
						flow = append(flow, distUpdate...)
					}
					if throwOut {
						break
					}
					if collision != nil {
						if m.Parameters.Volumetric {
							particlePtr.updateExtraDims(m)
							if particlePtr.y*particlePtr.y+particlePtr.z*particlePtr.z > m.Parameters.CathodeRadius*m.Parameters.CathodeRadius {
								break
							}
						}

						energyLoss := collision.Threshold
						cosChiScattered := m.Parameters.CrossSectionsData().SampleScatteringAngleCos(particlePtr.eKinetic, energyLoss, collision.Type)
						phi := 2. * math.Pi * rand.Float64()
						particlePtr.eKinetic -= energyLoss
						switch collision.Type {
						case lxgata.ELASTIC:

							// cosChiScattered = max(-1., min(1., (2.+particlePtr.eKinetic-2.*math.Pow(1.+particlePtr.eKinetic, rand.Float64()))/particlePtr.eKinetic))
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
							if m.Parameters.ParallelPlaneHollowCathode {
								const OMEGA = 15.
								ejected.eKinetic = OMEGA * math.Tan(rand.Float64()*math.Atan(particlePtr.eKinetic/(2.*OMEGA)))
							} else {
								ejected.eKinetic = particlePtr.eKinetic * rand.Float64()
							}
							cosChiEjected := math.Sqrt(ejected.eKinetic / particlePtr.eKinetic)
							ejected.redirect(cosChiEjected, math.Cos(phi+math.Pi), m)
							if ejected.totEnergy >= ionizationThreshold || m.Parameters.CalculateDistribution {
								computeWg.Add(1)
								computeflow <- &ejected
							}

							eScattered := particlePtr.eKinetic - ejected.eKinetic
							cosChiScattered = math.Sqrt(eScattered / particlePtr.eKinetic)
							particlePtr.eKinetic = eScattered
						}
						collflow <- CollisionEvent{int(particlePtr.x / m.XStep), energyLoss, collision.Type, particlePtr.origin}
						if collision.Type == lxgata.ATTACHMENT {
							break
						}
						particlePtr.redirect(cosChiScattered, math.Cos(phi), m)
					} else {
						currentCell := int(particlePtr.x / m.XStep)
						if currentCell < m.NumCells && m.Parameters.CountNulls {
							select {
							case collflow <- CollisionEvent{currentCell, 0, "NULL", particlePtr.origin}:
							default:
							}

						}
					}
				}
				if particlePtr.totEnergy < ionizationThreshold {
					ooeflow <- int(particlePtr.x / m.XStep)
				}
				if m.Parameters.CalculateDistribution {
					distflow <- flow
				}
				computeWg.Done()
			}
		}()
	}
	computeWg.Wait()
	close(computeflow)
	close(collflow)
	if m.Parameters.CalculateDistribution {
		close(distflow)
	}
	close(ooeflow)
	stateWg.Wait()
	print("\r")
}
