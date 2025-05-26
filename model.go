package main

import (
	"fmt"
	"math"
	"math/rand"
	"sync"

	"github.com/wildstyl3r/lxgata"
)

type Model struct {
	// crossSections *lxgata.Collisions
	parameters ModelParameters
	Vc         float64 // cathode fall potential
	// gasDensity               float64 // N
	inverseAbsConstEField    float64
	inverseNegativeVc        float64
	inverseCathodeFallLength float64

	Va float64 // additional voltage to avoid numeric negative energy beyond cathode fall

	numCells int

	xStep float64
	eStep float64

	lookUpPotential         []float64
	lookupTotalCrossSection []float64 // CS[energy]
	energyGrid              []float64
	gridBoundSqrt           []float64

	collisionAtCell     map[lxgata.CollisionType][][]uint16
	energyLossByProcess map[lxgata.CollisionType][]float64

	outOfEnergyAtCell []int
}

func newModel(CathodeFallLength float64, parameters ModelParameters) Model {
	m := Model{}
	m.parameters = parameters
	m.parameters.CathodeFallLength = CathodeFallLength
	m.Vc = math.Abs(m.parameters.CathodeFallPotential)
	m.Va = math.Abs(m.parameters.ConstEField * (m.parameters.GapLength - m.parameters.CathodeFallLength))
	m.inverseAbsConstEField = math.Abs(1. / m.parameters.ConstEField)
	m.inverseNegativeVc = -1. / m.Vc
	m.inverseCathodeFallLength = 1. / m.parameters.CathodeFallLength

	meanFreePath := 1. / (m.parameters._crossSections.SurplusCrossSection() * m.parameters.gasDensity)
	if m.parameters._verbose {
		fmt.Printf("Mean free path: %f\n", meanFreePath)
	}

	m.numCells = 5 * int(m.parameters.GapLength/meanFreePath)
	if m.parameters._verbose {
		fmt.Printf("xCells: %d;\n", m.numCells)
	}

	m.xStep = m.parameters.GapLength / float64(m.numCells+1)
	m.eStep = 0.01 // eV

	m.collisionAtCell = make(map[lxgata.CollisionType][][]uint16)
	m.energyLossByProcess = make(map[lxgata.CollisionType][]float64)
	for _, process := range []lxgata.CollisionType{
		lxgata.IONIZATION,
		lxgata.EXCITATION,
		lxgata.ELASTIC,
		lxgata.EFFECTIVE,
		lxgata.ATTACHMENT,
		lxgata.ROTATION,
		lxgata.CollisionType("NULL")} {
		if m.parameters.CalculateStdError {
			m.energyLossByProcess[process] = make([]float64, m.numCells+1)
			m.collisionAtCell[process] = make([][]uint16, m.numCells+1)
			for c := range m.numCells + 1 {
				m.collisionAtCell[process][c] = make([]uint16, parameters.NElectrons)
			}
		} else {
			m.energyLossByProcess[process] = make([]float64, m.numCells+1)
			m.collisionAtCell[process] = make([][]uint16, m.numCells+1)
			for c := range m.numCells + 1 {
				m.collisionAtCell[process][c] = make([]uint16, 1)
			}
		}
	}
	m.lookUpPotential = make([]float64, m.numCells+2)
	for cell := range m.lookUpPotential {
		m.lookUpPotential[cell] = m.VfromL(float64(cell) * m.xStep)
	}

	m.outOfEnergyAtCell = make([]int, m.numCells)

	numEnergyCells := int((m.Vc + m.Va + 5 + 0.1) / m.eStep)
	m.energyGrid = make([]float64, numEnergyCells)
	m.gridBoundSqrt = make([]float64, numEnergyCells)
	for i := range m.energyGrid {
		m.energyGrid[i] = float64(i) * m.eStep
		m.gridBoundSqrt[i] = math.Sqrt(m.energyGrid[i])
	}

	m.lookupTotalCrossSection = make([]float64, numEnergyCells)
	for i := range m.lookupTotalCrossSection {
		m.lookupTotalCrossSection[i] = m.parameters._crossSections.TotalCrossSectionAt(m.energyGrid[i])
	}
	return m
}

// g
func (s *Model) VfromL(l float64) (V float64) {
	cathodeFallPortion := l * s.inverseCathodeFallLength
	if l < s.parameters.CathodeFallLength {
		V = s.Vc * (cathodeFallPortion*(2.-cathodeFallPortion) - 1.)
	} else {
		V = -s.parameters.ConstEField * (l - s.parameters.CathodeFallLength)
	}
	V -= s.Va
	return V
}

// g^{-1}
func (s *Model) LfromV(V float64) (l float64) {
	if V < -s.Va {
		l = math.FMA(-s.parameters.CathodeFallLength, math.Sqrt((V+s.Va)*s.inverseNegativeVc), s.parameters.CathodeFallLength)
	} else {
		l = math.FMA((V + s.Va), s.inverseAbsConstEField, s.parameters.CathodeFallLength)
	}
	if math.IsNaN(l) {
		panic("x is NaN")
	}
	return
}

func (s *Model) EFieldFromL(l float64) (E float64) { // V/m
	E = s.parameters.ConstEField
	if l < s.parameters.CathodeFallLength {
		E = math.FMA(-2.*(1.-l*s.inverseCathodeFallLength)*s.Vc, s.inverseCathodeFallLength, E)
	}
	return
}

func (s *Model) EFieldFromPotential(V float64) (E float64) {
	E = s.parameters.ConstEField
	if V < -s.Va {
		E = math.FMA(-2*s.inverseCathodeFallLength, math.Sqrt(-(V+s.Va)*s.Vc), E)
	}
	return
}

func (m *Model) collisionSelector(eKinetic, totalCrossSectionPrimed float64) *lxgata.Collision {
	var crossSections = m.parameters._crossSections.CrossSectionsAt(eKinetic)
	var totalCrossSection = sumSlice(crossSections)
	var crossSectionAccum float64 = totalCrossSectionPrimed - totalCrossSection // the difference is null collision
	var choice float64 = rand.Float64() * totalCrossSectionPrimed

	selected := false
	if choice < crossSectionAccum {
		return nil
	}
	for i := range *m.parameters._crossSections {
		crossSectionAccum += crossSections[i]
		if choice < crossSectionAccum && !selected {
			return &(*m.parameters._crossSections)[i]
		}
	}
	return nil
}

func (m *Model) nextCollision(p *Particle) *lxgata.Collision {
	R := -math.Log(1. - rand.Float64())

	mt := m.parameters._crossSections.MinThresholdOfKind(lxgata.IONIZATION)
	energyCellIndex := int(p.eKinetic / m.eStep)
	maxEnergyCellIndex := int((m.Vc + m.Va + 5 + 0.1) / m.eStep)
	var CFLWorkaround bool = false
	var cachedVelocity float64
	var isVelCached = false
	// potential at dc = -Va
	for 0 <= energyCellIndex && energyCellIndex < maxEnergyCellIndex && p.totEnergy > mt {
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
		M *= m.parameters.gasDensity

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
					// p.x = m.LfromV(-(p.totEnergy - p.eStar))
					continue
				}
				if !CFLWorkaround {
					energyCellIndex--
				}
			} else {
				if eNext < p.eStar {
					println("k<*")
				}
				if !CFLWorkaround {
					energyCellIndex++
				}
			}
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
			if eColl > m.Vc+m.Va+5 || p.x < 0 || m.parameters.GapLength < p.x {
				break
			}
			return m.collisionSelector(eColl, totalCrossSectionPrimed)
		}
	}
	if p.totEnergy < mt {
		m.outOfEnergyAtCell[int(p.x/m.xStep)]++
	}
	p.x = m.xStep * float64(m.numCells+3)
	return nil
}

func (m *Model) run() {
	var computeWg, stateWg sync.WaitGroup

	collflow := make(chan CollisionEvent, 100000)
	if *m.parameters._dataFlags.all || *m.parameters._dataFlags.sequentials["Collision counters"].saveFlag || *m.parameters._dataFlags.sequentials["Normalized source term"].saveFlag || m.parameters.CalculateCathodeFallLength {
		stateWg.Add(1)
		go func() {
			for collision := range collflow {
				if collision.x < m.numCells+1 {
					m.collisionAtCell[collision.collType][collision.x][collision.origin]++
					if collision.collType != "NULL" {
						m.energyLossByProcess[collision.collType][collision.x] += collision.energyLoss
					}
				}
			}
			stateWg.Done()
		}()
	}

	computeflow := make(chan *Particle, m.parameters.NElectrons*100)
	for i := range m.parameters.NElectrons {
		origin := 0
		if m.parameters.CalculateStdError {
			origin = i
		}
		particle := m.newParticle(origin)
		computeWg.Add(1)
		computeflow <- &particle
	}

	status := []string{"//", "==", "\\\\", "||"}
	for range m.parameters._threads {
		go func() {
			counter := 0
			for particlePtr := range computeflow {
				counter++
				print("\r" + status[counter&0b11])
				mt := m.parameters._crossSections.MinThresholdOfKind(lxgata.IONIZATION)
				for (!m.parameters.ParallelPlaneHollowCathode && int((particlePtr.x+0.000001*m.xStep)/m.xStep) < m.numCells) ||
					(m.parameters.ParallelPlaneHollowCathode && int(particlePtr.x/m.xStep) < m.numCells+1) {
					if collision := m.nextCollision(particlePtr /*, stateflow*/); collision != nil {
						if m.parameters.Volumetric {
							particlePtr.updateExtraDims(m)
							if particlePtr.y*particlePtr.y+particlePtr.z*particlePtr.z > m.parameters.CathodeRadius*m.parameters.CathodeRadius {
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
							if m.parameters.ParallelPlaneHollowCathode {
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
						collflow <- CollisionEvent{int((particlePtr.x + m.xStep*0.5) / m.xStep), energyLoss, collision.Type, particlePtr.origin}
						if collision.Type == lxgata.ATTACHMENT {
							break
						}
						particlePtr.redirect(cosChiScattered, math.Cos(phi), m)
					} else {
						currentCell := int(particlePtr.x / m.xStep)
						if currentCell < m.numCells && m.parameters.CountNulls {
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
	stateWg.Wait()
	print("\r")
}
