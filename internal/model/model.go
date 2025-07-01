package model

import (
	"container/list"
	"fmt"
	"math"
	"math/rand"
	"sync"

	"github.com/wildstyl3r/lxgata"
	"github.com/wildstyl3r/stmc/internal/config"
	"github.com/wildstyl3r/stmc/internal/utils"
)

type Model struct {
	// crossSections *lxgata.Collisions
	Parameters config.ModelParameters
	Vc         float64 // cathode fall potential
	// gasDensity               float64 // N
	inverseAbsConstEField    float64
	inverseNegativeVc        float64
	inverseCathodeFallLength float64

	Va float64 // additional voltage to avoid numeric negative energy beyond cathode fall

	NumCells, NumMuCells int

	EStep, XStep, MuStep float64

	lookUpPotential         []float64
	lookupTotalCrossSection []float64 // CS[energy]
	LookupEnergy            []float64
	LookupInverseVelocity   []float64
	GridInverseMu           []float64
	lookupMNumerator        []float64

	Distribution [][]*list.List

	CollisionAtCell     map[lxgata.CollisionType][][]uint16
	EnergyLossByProcess map[lxgata.CollisionType][]float64

	OutOfEnergyAtCell []int
}

func NewModel(CathodeFallLength float64, parameters config.ModelParameters) Model {
	m := Model{}
	m.Parameters = parameters
	m.Parameters.CathodeFallLength = CathodeFallLength
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
	if m.Parameters.Verbose() {
		fmt.Printf("xCells: %d;\n", m.NumCells)
	}

	m.XStep = m.Parameters.GapLength / float64(m.NumCells+1)
	m.EStep = m.Parameters.EnergyStep // eV
	m.MuStep = m.Parameters.MuStep    // eV
	m.NumMuCells = int(2. / m.MuStep)

	m.CollisionAtCell = make(map[lxgata.CollisionType][][]uint16)
	m.EnergyLossByProcess = make(map[lxgata.CollisionType][]float64)
	for _, process := range []lxgata.CollisionType{
		lxgata.IONIZATION,
		lxgata.EXCITATION,
		lxgata.ELASTIC,
		lxgata.EFFECTIVE,
		lxgata.ATTACHMENT,
		lxgata.ROTATION} {
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

	numEnergyCells := int((m.Vc + m.Va + 5 + 0.1) / m.EStep)
	m.lookupTotalCrossSection = make([]float64, numEnergyCells)
	m.LookupEnergy = make([]float64, numEnergyCells)
	m.LookupInverseVelocity = make([]float64, numEnergyCells)
	m.lookupMNumerator = make([]float64, numEnergyCells)
	for gridNode := range m.LookupEnergy {
		m.LookupEnergy[gridNode] = float64(gridNode) * m.EStep
		m.lookupTotalCrossSection[gridNode] = m.Parameters.CrossSectionsData().TotalCrossSectionAt(m.LookupEnergy[gridNode])
		m.LookupInverseVelocity[gridNode] = 1. / utils.EV2electronVelocity(m.LookupEnergy[gridNode]+0.5*parameters.EnergyStep)
		m.lookupMNumerator[gridNode] = m.Parameters.GasDensity * m.lookupTotalCrossSection[gridNode] * math.Sqrt(m.LookupEnergy[gridNode])
	}
	m.GridInverseMu = make([]float64, int(1./m.MuStep)+1)
	for muNode := range m.GridInverseMu {
		m.GridInverseMu[muNode] = 1. / (m.MuStep * (float64(muNode) + 0.5))
	}

	if m.Parameters.CalculateDistribution {
		m.Distribution = make([][]*list.List, m.NumCells)
		for x := range m.Distribution {
			m.Distribution[x] = make([]*list.List, m.NumMuCells) //list.New()
			for mu := range m.NumMuCells {
				m.Distribution[x][mu] = list.New()
			}
		}
	}

	return m
}

func (m *Model) collisionSelector(eKinetic, x, M float64) *lxgata.Collision {
	var crossSections = m.Parameters.CrossSectionsData().CrossSectionsAt(eKinetic)
	var totalCrossSection, totalCrossSectionPrimed = utils.SumSlice(crossSections), M * math.Abs(m.EFieldFromL(x)) / math.Sqrt(eKinetic) / m.Parameters.GasDensity
	var crossSectionAccum float64 = totalCrossSectionPrimed - totalCrossSection // the difference is null collision
	var choice float64 = rand.Float64() * totalCrossSectionPrimed

	selected := false
	if choice < crossSectionAccum {
		return nil
	}
	for i := range *m.Parameters.CrossSectionsData() {
		crossSectionAccum += crossSections[i]
		if choice < crossSectionAccum && !selected {
			return &(*m.Parameters.CrossSectionsData())[i]
		}
	}
	return nil
}

type FlowElem struct {
	x, mu  int
	energy float64
}

type Flow []FlowElem

func (m *Model) getFlowBetweenEnergies(startEnergy, endEnergy float64, p *Particle) (flow []FlowElem) {
	energyMin, energyMax := min(startEnergy, endEnergy), max(startEnergy, endEnergy)
	var direction float64
	if startEnergy < endEnergy {
		direction = +1.
	} else {
		direction = -1.
	}
	xMin, xMax := m.LfromV(-(p.totEnergy - energyMin)), m.LfromV(-(p.totEnergy - energyMax))
	gridMinIndex, gridMaxIndex := int(xMin/m.XStep)+1, int(xMax/m.XStep)+1
	if gridMinIndex < 0 {
		panic("i<0")
	}
	for i := gridMinIndex; i < gridMaxIndex; i++ {
		eKinetic := p.totEnergy + m.VfromL(float64(i)*m.XStep)
		mu := math.Copysign(math.Sqrt((eKinetic-p.eStar)/eKinetic), direction)
		flow = append(flow, FlowElem{i, int(math.Floor(mu/m.MuStep)) + m.NumMuCells/2, eKinetic})
	}
	return
}

func (m *Model) nextCollision(p *Particle) (collisionType *lxgata.Collision, flow []FlowElem, throwOut bool) {
	R := -math.Log(1. - rand.Float64())

	minEnergy, maxEnergy := max(0, p.totEnergy-(m.Vc+m.Va)), p.totEnergy
	var currentCellIndex int
	minEnergyCellIndex, maxEnergyCellIndex := int(minEnergy/m.EStep), int(maxEnergy/m.EStep)+1
	// sheathLengthEnergy := p.totEnergy - m.VfromL(m.Parameters.CathodeFallLength)
	// sheathLengthEnergyCell := int(sheathLengthEnergy/m.EStep)
	// var CFLWorkaround bool = false
	// var cachedVelocity float64
	// var isVelCached = false
	// potential at dc = -Va
	// flow = make([]FlowElem, 0)

	alignedToEnergyGrid := false

	getPossibleCollision := func(segmentStartEnergy, M float64) (collisionType *lxgata.Collision, reversalHappened bool, throwOut bool) {
		var delta float64
		if p.mu < 0 {
			delta = math.Sqrt(segmentStartEnergy-p.eStar) - R/(2.*M) // == sqrt(pColl - p.eStar) i.e. abs(speed) equivalent along x axis
			if delta < 0 {
				R -= 2 * M * math.Sqrt(segmentStartEnergy-p.eStar)
				p.mu = +0.
				p.eKinetic = p.eStar
				reversalHappened = true
				return
			}
		} else {
			delta = R/(2.*M) + math.Sqrt(segmentStartEnergy-p.eStar)
		}
		collisionEnergy := math.FMA(delta, delta, p.eStar) // R = 2M[sqrt(p.e-p.eStar) - sqrt(eColl - p.eStar)]

		if m.Parameters.CalculateDistribution {
			flow = append(flow, m.getFlowBetweenEnergies(segmentStartEnergy, collisionEnergy, p)...)
		}
		p.setEnergy(collisionEnergy, m, true, true)
		if maxEnergy < collisionEnergy || p.x < 0 || m.Parameters.GapLength < p.x {
			throwOut = true
			return
		}
		return m.collisionSelector(collisionEnergy, p.x, M), false, false
	}

	// aroundSheathLength := false

	for !alignedToEnergyGrid || (minEnergyCellIndex <= currentCellIndex && currentCellIndex <= maxEnergyCellIndex && 0 <= p.x && p.x < m.Parameters.GapLength) {
		var lowEnergy, segmentStartEnergy, segmentEndEnergy, highEnergy float64
		var nextCellIndex, highEnergyCellIndex int
		var reversalHappened, arrivalAtCathode, arrivalAtGapEnd, highEnergyAligned bool
		if p.mu < 0 {
			if !alignedToEnergyGrid {
				nextCellIndex = int(p.eKinetic / m.EStep)
				lowEnergy, highEnergy = m.LookupEnergy[nextCellIndex], p.eKinetic

			} else {
				if currentCellIndex <= minEnergyCellIndex {
					throwOut = true
					return
				}
				lowEnergy, highEnergy = m.LookupEnergy[currentCellIndex-1], m.LookupEnergy[currentCellIndex]
				highEnergyCellIndex, highEnergyAligned = currentCellIndex, true
				nextCellIndex = currentCellIndex - 1
			}

			//check for reversal
			if lowEnergy < p.eStar {
				lowEnergy, alignedToEnergyGrid = p.eStar, false
				nextCellIndex = currentCellIndex
				reversalHappened = true
			}

			//check for cathode return
			if p.totEnergy-lowEnergy > m.Vc+m.Va {
				arrivalAtCathode = true
				lowEnergy = p.totEnergy - (m.Vc + m.Va)
			}
			segmentStartEnergy, segmentEndEnergy = highEnergy, lowEnergy
		} else {
			if !alignedToEnergyGrid {
				nextCellIndex = int(p.eKinetic/m.EStep) + 1
				lowEnergy, highEnergy = p.eKinetic, m.LookupEnergy[nextCellIndex]
			} else {
				if currentCellIndex >= maxEnergyCellIndex {
					throwOut = true
					return
				}
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
		G = 2 * M * (math.Sqrt(highEnergy-p.eStar) - math.Sqrt(lowEnergy-p.eStar))

		if G < R {
			if m.Parameters.CalculateDistribution {
				flow = append(flow, m.getFlowBetweenEnergies(segmentStartEnergy, segmentEndEnergy, p)...)
			}
			R -= G
		} else {
			var strangeReversalHappened bool
			collisionType, strangeReversalHappened, throwOut = getPossibleCollision(segmentStartEnergy, M)
			if !strangeReversalHappened {
				return
			} else {
				reversalHappened = true
			}
		}

		if arrivalAtCathode || arrivalAtGapEnd {
			throwOut = true
			return
		}

		if reversalHappened {
			p.mu = +0.
			p.eKinetic = p.eStar
			alignedToEnergyGrid = false
		} else {
			alignedToEnergyGrid = true
		}
		currentCellIndex = nextCellIndex

		// if cellIndex == m.numCells && m.parameters.ParallelPlaneHollowCathode {
		// 	cellIndex--
		// 	p.mu = -p.mu
		// }

		// if (p.totEnergy-leftBound) > m.Va &&
		// 	m.Va > (p.totEnergy-rightBound) &&
		// 	!CFLWorkaround {
		// 	CFLWorkaround = true
		// 	// isVelCached = false
		// 	if (p.totEnergy - p.eKinetic) > m.Va {
		// 		rightBound = p.totEnergy - m.Va
		// 	} else {
		// 		leftBound = p.totEnergy - m.Va
		// 	}
		// } else {
		// 	CFLWorkaround = false
		// }

	}
	// if p.totEnergy < ionizationThreshold && p.x < m.Parameters.GapLength {
	// 	m.OutOfEnergyAtCell[int(p.x/m.XStep)]++
	// }
	return nil, flow, false
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
					if distUpdate[i].x < m.NumCells {
						m.Distribution[distUpdate[i].x][distUpdate[i].mu].PushBack(distUpdate[i].energy)
					}
				}
			}
			stateWg.Done()
		}()
	}

	var ooeflow chan int
	stateWg.Add(1)
	go func() {
		for ooeEvent := range ooeflow {
			if ooeEvent < m.NumCells {
				m.OutOfEnergyAtCell[ooeEvent]++
			}
		}
		stateWg.Done()
	}()

	computeflow := make(chan *Particle, m.Parameters.NElectrons*100)
	for i := range m.Parameters.NElectrons {
		origin := 0
		if m.Parameters.CalculateStdError {
			origin = i
		}
		particle := m.newParticle(origin)
		if m.Parameters.CalculateDistribution {
			m.Distribution[0][int(math.Floor(particle.mu/m.MuStep))+m.NumMuCells/2].PushBack(particle.eKinetic)
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
				for ionizationThreshold < particlePtr.totEnergy || m.Parameters.CalculateDistribution { //xCell :=int((particlePtr.x+0.5*m.XStep)/m.XStep); (!m.Parameters.ParallelPlaneHollowCathode && xCell < m.NumCells) ||
					// (m.Parameters.ParallelPlaneHollowCathode && int(particlePtr.x/m.XStep) < m.NumCells+1) {
					collision, distUpdate, throwOut := m.nextCollision(particlePtr /*, stateflow*/)
					if m.Parameters.CalculateDistribution {
						distflow <- distUpdate
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
						cosChiScattered := 1. - 2.*rand.Float64()
						phi := 2. * math.Pi * rand.Float64()
						energyLoss := collision.Threshold
						particlePtr.eKinetic -= energyLoss
						switch collision.Type {
						case lxgata.ELASTIC:
							cosChiScattered = max(-1., min(1., (2.+particlePtr.eKinetic-2.*math.Pow(1.+particlePtr.eKinetic, rand.Float64()))/particlePtr.eKinetic))
							energyLoss = particlePtr.eKinetic * 2. * collision.MassRatio * (1. - cosChiScattered)
							particlePtr.eKinetic -= energyLoss

						case lxgata.EFFECTIVE:
							energyLoss = particlePtr.eKinetic * 2. * collision.MassRatio * (1. - cosChiScattered)
							particlePtr.eKinetic -= energyLoss

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
