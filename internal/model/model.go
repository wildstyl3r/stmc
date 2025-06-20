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
	// crossSections *lxgata.Collisions
	Parameters config.ModelParameters
	Vc         float64 // cathode fall potential
	// gasDensity               float64 // N
	inverseAbsConstEField    float64
	inverseNegativeVc        float64
	inverseCathodeFallLength float64

	Va float64 // additional voltage to avoid numeric negative energy beyond cathode fall

	NumCells int

	XStep float64
	EStep float64

	lookUpPotential         []float64
	lookupTotalCrossSection []float64 // CS[energy]
	energyGrid              []float64
	gridBoundSqrt           []float64

	Distribution [][]int

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
	m.energyGrid = make([]float64, numEnergyCells)
	m.gridBoundSqrt = make([]float64, numEnergyCells)
	for gridNode := range m.energyGrid {
		m.energyGrid[gridNode] = float64(gridNode) * m.EStep
		m.gridBoundSqrt[gridNode] = math.Sqrt(m.energyGrid[gridNode])
	}

	m.lookupTotalCrossSection = make([]float64, numEnergyCells)
	for i := range m.lookupTotalCrossSection {
		m.lookupTotalCrossSection[i] = m.Parameters.CrossSectionsData().TotalCrossSectionAt(m.energyGrid[i])
	}

	m.Distribution = make([][]int, int((m.Vc+m.Va+10)/parameters.EnergyStep))
	for i := range m.Distribution {
		m.Distribution[i] = make([]int, m.NumCells)
	}
	return m
}

func (m *Model) collisionSelector(eKinetic, totalCrossSectionPrimed float64) *lxgata.Collision {
	var crossSections = m.Parameters.CrossSectionsData().CrossSectionsAt(eKinetic)
	var totalCrossSection = utils.SumSlice(crossSections)
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
	xStart, xEnd, e, sign int
}

type Flow struct {
	seq []FlowElem
}

func (m *Model) nextCollision(p *Particle) (*lxgata.Collision, Flow) {
	R := -math.Log(1. - rand.Float64())

	mt := m.Parameters.CrossSectionsData().MinThresholdOfKind(lxgata.IONIZATION)
	energyCellIndex := int(p.eKinetic / m.EStep)
	alignedToEnergyGrid := false
	maxEnergyCellIndex := int((m.Vc + m.Va + 5 + 0.1) / m.EStep)
	var CFLWorkaround bool = false
	var cachedVelocity float64
	var isVelCached = false
	// potential at dc = -Va
	flow := Flow{seq: make([]FlowElem, 0)}
	for 0 <= energyCellIndex && energyCellIndex < maxEnergyCellIndex && (p.totEnergy > mt || true) && p.x < m.Parameters.GapLength {
		// if cellIndex == m.numCells && m.parameters.ParallelPlaneHollowCathode {
		// 	cellIndex--
		// 	p.mu = -p.mu
		// }

		leftBound, rightBound := m.energyGrid[energyCellIndex], m.energyGrid[energyCellIndex+1]

		var M, totalCrossSectionPrimed float64
		if (p.totEnergy-leftBound) > m.Va &&
			m.Va > (p.totEnergy-rightBound) &&
			!CFLWorkaround {
			CFLWorkaround = true
			isVelCached = false
			if (p.totEnergy - p.eKinetic) > m.Va {
				rightBound = p.totEnergy - m.Va
			} else {
				leftBound = p.totEnergy - m.Va
			}
		} else {
			CFLWorkaround = false
		}

		M = m.lookupTotalCrossSection[energyCellIndex+1] * m.gridBoundSqrt[energyCellIndex+1] / -m.EFieldFromPotential(rightBound-p.totEnergy)
		totalCrossSectionPrimed = -M * m.EFieldFromPotential(p.eKinetic-p.totEnergy) / math.Sqrt(p.eKinetic)
		M *= m.Parameters.GasDensity

		var G, eNext float64
		if p.mu < 0 {
			tempVel := math.Sqrt(leftBound - p.eStar)
			if isVelCached {
				G = 2 * M * (cachedVelocity - tempVel)
			} else {
				G = 2 * M * (math.Sqrt(p.eKinetic-p.eStar) - tempVel)
			}
			cachedVelocity = tempVel
			isVelCached = true
			// G = 2 * M * (math.Sqrt(p.eKinetic-p.eStar) - math.Sqrt(leftBound-p.eStar))
			eNext = leftBound
		} else {
			tempVel := math.Sqrt(rightBound - p.eStar)
			if isVelCached {
				G = 2 * M * (tempVel - cachedVelocity)
			} else {
				G = 2 * M * (tempVel - math.Sqrt(p.eKinetic-p.eStar))
			}
			cachedVelocity = tempVel
			isVelCached = true
			// G = 2 * M * (math.Sqrt(rightBound-p.eStar) - math.Sqrt(p.eKinetic-p.eStar))
			eNext = rightBound
		}

		if G < R { // no collision
			if p.mu < 0 {
				if eNext < p.eStar {
					R -= 2 * M * math.Sqrt(p.eKinetic-p.eStar)
					p.mu = +0.
					p.eKinetic = p.eStar
					isVelCached = false
					alignedToEnergyGrid = false
					// p.x = m.LfromV(-(p.totEnergy - p.eStar))
					continue
				}
				if !CFLWorkaround {
					if alignedToEnergyGrid {
						flow.seq = append(flow.seq, FlowElem{xStart: int(m.LfromV(-(p.totEnergy - (float64(energyCellIndex)+0.0)*m.Parameters.EnergyStep)) / m.XStep), xEnd: int(m.LfromV(-(p.totEnergy - (float64(energyCellIndex-1)+0.0)*m.Parameters.EnergyStep)) / m.XStep), e: int(((float64(energyCellIndex)+0.5)*m.Parameters.EnergyStep-p.eStar)/m.Parameters.EnergyStep) - 1, sign: -1})
					} else {
						flow.seq = append(flow.seq, FlowElem{xStart: int(p.x / m.XStep), xEnd: int(m.LfromV(-(p.totEnergy - (float64(energyCellIndex-1)+0.0)*m.Parameters.EnergyStep)) / m.XStep), e: int(((float64(energyCellIndex)+0.5)*m.Parameters.EnergyStep-p.eStar)/m.Parameters.EnergyStep) - 1, sign: -1})
					}

					energyCellIndex--
				}
			} else {
				if eNext < p.eStar {
					println("k<*")
				}
				if !CFLWorkaround {
					energyCellIndex++
					if alignedToEnergyGrid {
						flow.seq = append(flow.seq, FlowElem{xStart: int(m.LfromV(-(p.totEnergy - (float64(energyCellIndex-1)+0.0)*m.Parameters.EnergyStep)) / m.XStep), xEnd: int(m.LfromV(-(p.totEnergy - (float64(energyCellIndex)+0.0)*m.Parameters.EnergyStep)) / m.XStep), e: int(((float64(energyCellIndex)+0.5)*m.Parameters.EnergyStep-p.eStar)/m.Parameters.EnergyStep) - 1, sign: +1})
					} else {
						flow.seq = append(flow.seq, FlowElem{xStart: int(m.LfromV(-(p.totEnergy - (float64(energyCellIndex-1)+0.0)*m.Parameters.EnergyStep)) / m.XStep), xEnd: int(p.x / m.XStep), e: int(((float64(energyCellIndex)+0.5)*m.Parameters.EnergyStep-p.eStar)/m.Parameters.EnergyStep) - 1, sign: +1})
					}

				}
			}
			alignedToEnergyGrid = true
			if eNext > m.Vc+m.Va+5 {
				break
			}
			p.setEnergy(eNext, m, true, false)
			R -= G
		} else { // collision occured (possibly null)
			var eColl float64
			if p.mu < 0 {
				delta := math.Sqrt(p.eKinetic-p.eStar) - R/(2.*M) // == sqrt(pColl - p.eStar) i.e. abs(speed) equivalent along x axis
				if delta >= 0 {
					eColl = math.FMA(delta, delta, p.eStar) // R = 2M[sqrt(p.e-p.eStar) - sqrt(eColl - p.eStar)]
				} else {
					R -= 2. * M * (math.Sqrt(p.eKinetic - p.eStar))
					p.mu = +0.
					p.eKinetic = p.eStar
					isVelCached = false
					// p.x = m.LfromV(-(p.totEnergy - p.eStar))
					continue
				}
			} else {
				t := R/(2.*M) + math.Sqrt(p.eKinetic-p.eStar)
				eColl = math.FMA(t, t, p.eStar)
			}
			p.setEnergy(eColl, m, true, true)
			if eColl > m.Vc+m.Va+5 || p.x < 0 || m.Parameters.GapLength < p.x {
				break
			}
			return m.collisionSelector(eColl, totalCrossSectionPrimed), flow
		}
	}
	if p.totEnergy < mt && p.x < m.Parameters.GapLength {
		m.OutOfEnergyAtCell[int(p.x/m.XStep)]++
	}
	p.x = m.XStep * float64(m.NumCells+3)
	return nil, flow
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

	distflow := make(chan Flow, 10000)
	stateWg.Add(1)
	go func() {
		for distUpdate := range distflow {
			for i := range distUpdate.seq {
				for x := max(distUpdate.seq[i].xStart, 0); x < distUpdate.seq[i].xEnd && x < m.NumCells; x++ {
					if 0 <= distUpdate.seq[i].e && distUpdate.seq[i].e < len(m.Distribution) {
						m.Distribution[distUpdate.seq[i].e][x] += distUpdate.seq[i].sign
					}
				}
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
		m.Distribution[int(particle.eKinetic/m.Parameters.EnergyStep)][0]++
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
				mt := m.Parameters.CrossSectionsData().MinThresholdOfKind(lxgata.IONIZATION)
				for (!m.Parameters.ParallelPlaneHollowCathode && int((particlePtr.x+0.000001*m.XStep)/m.XStep) < m.NumCells) ||
					(m.Parameters.ParallelPlaneHollowCathode && int(particlePtr.x/m.XStep) < m.NumCells+1) {
					collision, distUpdate := m.nextCollision(particlePtr /*, stateflow*/)
					distflow <- distUpdate
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
							cosChiScattered = (2. + particlePtr.eKinetic - 2.*math.Pow(1.+particlePtr.eKinetic, rand.Float64())) / particlePtr.eKinetic
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
							if ejected.totEnergy >= mt {
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
				computeWg.Done()
			}
		}()
	}
	computeWg.Wait()
	close(computeflow)
	close(collflow)
	close(distflow)
	stateWg.Wait()
	print("\r")
}
