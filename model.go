package main

import (
	"fmt"
	"math"
	"math/rand"
	"sync"

	"github.com/wildstyl3r/lxgata"
)

type Model struct {
	crossSections lxgata.Collisions
	parameters    ModelParameters
	dataFlags     DataFlags
	Vc            float64 // cathode fall potential
	gasDensity    float64 // N
	threads       int

	Va float64 // additional voltage to avoid numeric negative energy beyond cathode fall

	numCells   int
	numCellsE  int
	numCellsMu int

	xStep float64

	lookUpPotential []float64

	distribution [][][]int // psi(x,e,mu)

	collisionAtCell     map[string][]int
	energyLossByProcess map[string][]float64

	probabilities map[string][]AveragingElement

	// electronsAtCell []int
}

func (m *Model) init() {
	m.Va = -m.parameters.ConstEField * (m.parameters.GapLength - m.parameters.CathodeFallLength)

	meanFreePath := 1. / (m.crossSections.SurplusCrossSection() * m.gasDensity)
	fmt.Printf("Mean free path: %f\n", meanFreePath)

	m.numCells = 5 * int(m.parameters.GapLength/meanFreePath)
	m.numCellsE = int((math.Abs(m.parameters.CathodeFallPotential)+math.Abs(m.Va)+1)/m.parameters.DeltaE) + 1
	m.numCellsMu = int(2./m.parameters.DeltaMu) + 1
	fmt.Printf("xCells: %d; eCells: %d; muCells: %d\n", m.numCells, m.numCellsE, m.numCellsMu)

	m.xStep = m.parameters.GapLength / float64(m.numCells+1)

	m.distribution = make([][][]int, m.numCells+1)
	m.collisionAtCell = make(map[string][]int)
	m.energyLossByProcess = make(map[string][]float64)
	m.probabilities = make(map[string][]AveragingElement)
	// m.electronsAtCell = make([]int, m.numCells+1)
	m.lookUpPotential = make([]float64, m.numCells+2)

	for cell := range m.distribution {
		m.distribution[cell] = make([][]int, m.numCellsE)
		for energyCell := range m.distribution[cell] {
			m.distribution[cell][energyCell] = make([]int, m.numCellsMu)
		}
	}
	for cell := range m.lookUpPotential {
		m.lookUpPotential[cell] = m.VfromL(float64(cell) * m.xStep)
	}
}

// g
func (s *Model) VfromL(l float64) (V float64) {
	cathodeFallPortion := l / s.parameters.CathodeFallLength
	if l < s.parameters.CathodeFallLength {
		V = math.Abs(s.Vc) * (cathodeFallPortion*(2.-cathodeFallPortion) - 1.)
	} else {
		V = -s.parameters.ConstEField * (l - s.parameters.CathodeFallLength)
	}
	V -= math.Abs(s.Va)
	return V
}

// g^-1
func (s *Model) LfromV(V float64) (l float64) {
	Vtot := math.Abs(s.Va) + math.Abs(s.Vc)
	if V < -math.Abs(s.Va) {
		l = s.parameters.CathodeFallLength * (1. - math.Sqrt(1.-(V+Vtot)/math.Abs(s.Vc)))
	} else {
		l = (V+math.Abs(s.Va))/(-s.parameters.ConstEField) + s.parameters.CathodeFallLength
	}
	if math.IsNaN(l) {
		panic("x is NaN")
	}
	return l
}

func (s *Model) EFieldFromL(l float64) (E float64) { // V/m
	E = s.parameters.ConstEField
	if l < s.parameters.CathodeFallLength {
		E += -2. * (1. - l/s.parameters.CathodeFallLength) * math.Abs(s.Vc) / s.parameters.CathodeFallLength
	}
	return E
}

func (s *Model) collisionSelector(p *Particle, eKinetic, M float64, probChan chan<- ProbIncrement) *lxgata.Collision {
	var totalCrossSectionPrimed float64 = M * math.Abs(s.EFieldFromL(p.x)) / (s.gasDensity * math.Sqrt(eKinetic)) // [m^2]
	var crossSectionAccum float64 = totalCrossSectionPrimed - s.crossSections.TotalCrossSectionAt(eKinetic)       // the difference is null collision
	var choice float64 = rand.Float64()
	selector := crossSectionAccum / totalCrossSectionPrimed

	currentCell := int(p.x / s.xStep)
	eAtCellBound := s.lookUpPotential[currentCell] + p.totEnergy
	var totalCrossSectionAtCellBound float64 = s.crossSections.TotalCrossSectionAt(eAtCellBound)
	var totalCrossSectionPrimedAtCellBound float64 = M * math.Abs(s.EFieldFromL(p.x)) / (s.gasDensity * math.Sqrt(eAtCellBound)) // [m^2]
	probChan <- ProbIncrement{currentCell, "NULL", (totalCrossSectionPrimedAtCellBound - totalCrossSectionAtCellBound) / totalCrossSectionAtCellBound}

	var selection *lxgata.Collision
	selected := false
	if choice < selector {
		selection = nil // null collision
		selected = true
	}
	for i, collision := range s.crossSections {
		crossSectionAccum += collision.CrossSectionAt(eKinetic)
		selector = crossSectionAccum / totalCrossSectionPrimed
		probChan <- ProbIncrement{currentCell, string(collision.Type), collision.CrossSectionAt(eAtCellBound) / totalCrossSectionAtCellBound}
		if choice < selector && !selected {
			selection = &s.crossSections[i]
			selected = true
		}
	}
	return selection
}

func (s *Model) nextCollision(p *Particle, stateChan chan<- StateIncrement, probChan chan<- ProbIncrement) *lxgata.Collision { //ch chan<- StateIncrement,, acc map[StateIncrement]int
	R := -math.Log(1. - rand.Float64())

	cellIndex := int((p.x + 0.000001*s.xStep) / s.xStep)
	for (cellIndex < s.numCells && !s.parameters.ParallelPlaneHollowCathode) || (cellIndex < s.numCells+1 && s.parameters.ParallelPlaneHollowCathode) {
		if p.x < 0 || (p.x == 0 && p.mu < 0) || (s.parameters.ParallelPlaneHollowCathode && cellIndex == s.numCells && p.totEnergy < s.crossSections.MinThresholdOfKind(lxgata.IONIZATION)) { // must do this only after some time dut to density averaging
			p.x = s.xStep * float64(s.numCells+3)
			return nil
		}
		if cellIndex == s.numCells {
			cellIndex--
			p.mu = -p.mu
		}
		eLeft, eRight := max(p.totEnergy+s.lookUpPotential[cellIndex], p.eStar), p.totEnergy+s.lookUpPotential[cellIndex+1]
		if p.mu < 0 && p.eKinetic == eLeft && eLeft != p.eStar {
			cellIndex--
			eLeft, eRight = max(p.totEnergy+s.lookUpPotential[cellIndex], p.eStar), p.totEnergy+s.lookUpPotential[cellIndex+1]
		}
		M := p.M(cellIndex, s)
		var G, eNext float64
		if p.mu < 0 {
			G = 2. * M * (math.Sqrt(p.eKinetic-p.eStar) - math.Sqrt(eLeft-p.eStar))
			eNext = eLeft
		} else {
			G = 2. * M * (math.Sqrt(eRight-p.eStar) - math.Sqrt(p.eKinetic-p.eStar))
			eNext = eRight
		}
		if G < R { // no collision
			if p.eStar == eNext {
				p.mu = +0.
			}
			p.setEnergy(eNext, s, p.eStar == eNext)
			cellIndex = int((p.x + 0.000001*s.xStep) / s.xStep)
			{
				xIndex := int((p.x + s.xStep/2.) / s.xStep)
				eIndex := int(p.eKinetic / s.parameters.DeltaE)
				muIndex := int((p.mu + 1. + s.parameters.DeltaMu/2.) / s.parameters.DeltaMu)
				if eIndex < s.numCellsE && muIndex < s.numCellsMu && xIndex < s.numCells {
					stateChan <- StateIncrement{xIndex, eIndex, muIndex}
				}
			}
			R -= G
		} else { // collision occured (possibly null)
			var eColl float64
			if p.mu < 0 {
				delta := math.Sqrt(p.eKinetic-p.eStar) - R/(2.*M) // == sqrt(pColl - p.eStar) i.e. abs(speed) equivalent along x axis
				if delta >= 0 {
					eColl = math.Pow(delta, 2.) + p.eStar // R = 2M[sqrt(p.e-p.eStar) - sqrt(eColl - p.eStar)]
				} else {
					R -= 2. * M * (math.Sqrt(p.eKinetic - p.eStar))
					p.mu = +0.
					p.setEnergy(p.eStar, s, true) // passing through stop point
					cellIndex = int((p.x + 0.000001*s.xStep) / s.xStep)
					continue
				}
			} else {
				eColl = math.Pow(R/(2.*M)+math.Sqrt(p.eKinetic-p.eStar), 2.) + p.eStar // R = 2M[sqrt(eColl - p.eStar) - sqrt(p.e-p.eStar)]
			}
			p.setEnergy(eColl, s, true)
			return s.collisionSelector(p, eColl, M, probChan)
		}
	}
	return nil
}

func (m *Model) run() {
	var computeWg, stateWg sync.WaitGroup
	stateflow := make(chan StateIncrement, 100000) //make(chan map[StateIncrement]int, 1000)

	stateWg.Add(1)
	go func() {
		for update := range stateflow {
			value := 1

			if update.muIndex < m.numCellsMu/2 {
				m.distribution[update.xIndex][update.eIndex][update.muIndex] += value
			} else {
				m.distribution[update.xIndex][update.eIndex][update.muIndex] += value
			}
		}
		stateWg.Done()
	}()

	probflow := make(chan ProbIncrement, 100000)
	if *m.dataFlags.all || *m.dataFlags.specification["Actual process probabilities"].saveFlag {
		stateWg.Add(1)
		go func() {
			for probInc := range probflow {
				if m.probabilities[probInc.collType] == nil {
					m.probabilities[probInc.collType] = make([]AveragingElement, m.numCells+1)
				}
				m.probabilities[probInc.collType][probInc.x].add(probInc.prob)
			}
			stateWg.Done()
		}()
	}

	collflow := make(chan CollisionEvent, 100000)
	if *m.dataFlags.all || *m.dataFlags.specification["Collision counters"].saveFlag {
		stateWg.Add(1)
		go func() {
			for collision := range collflow {
				if m.collisionAtCell[collision.collType] == nil {
					m.collisionAtCell[collision.collType] = make([]int, m.numCells+1)
				}
				if m.energyLossByProcess[collision.collType] == nil {
					m.energyLossByProcess[collision.collType] = make([]float64, m.numCells+1)
				}
				m.collisionAtCell[collision.collType][collision.x]++
				m.energyLossByProcess[collision.collType][collision.x] += collision.energyLoss
			}
			stateWg.Done()
		}()
	}

	computeflow := make(chan *Particle, m.parameters.NElectrons*100)
	for i := 0; i < m.parameters.NElectrons; i++ {
		particle := m.newParticle()
		m.distribution[0][int(particle.eKinetic/m.parameters.DeltaE)][int((particle.mu+1.+m.parameters.DeltaMu/2.)/m.parameters.DeltaMu)]++
		computeflow <- &particle
		computeWg.Add(1)
	}

	status := []string{"//", "==", "\\\\", "||"}
	for w := 0; w < m.threads; w++ {
		go func() {
			counter := 0
			for particlePtr := range computeflow {
				counter++
				print("\r" + status[counter&0b11])
				for (!m.parameters.ParallelPlaneHollowCathode && int((particlePtr.x+0.000001*m.xStep)/m.xStep) < m.numCells) || (m.parameters.ParallelPlaneHollowCathode && int(particlePtr.x/m.xStep) < m.numCells+1) {
					if collision := m.nextCollision(particlePtr, stateflow, probflow); collision != nil {
						// fmt.Printf("x:%f; cell: %d of %d; mu: %f; eKin: %f\n", particlePtr.x, int(particlePtr.x/m.xStep), m.numCells, particlePtr.mu, particlePtr.eKinetic)
						cosChi := 1. - 2.*rand.Float64()
						phi := 2. * math.Pi * rand.Float64()
						switch collision.Type {
						case lxgata.ELASTIC, lxgata.EFFECTIVE:
							collflow <- CollisionEvent{int((particlePtr.x + m.xStep/2) / m.xStep), particlePtr.eKinetic * (1. - (2.*collision.MassRatio)*(1.-cosChi)), string(collision.Type)}
							particlePtr.eKinetic *= (1. - 2.*collision.MassRatio*(1.-cosChi))
							particlePtr.redirect(cosChi, math.Cos(phi), m)

						case lxgata.EXCITATION:
							collflow <- CollisionEvent{int((particlePtr.x + m.xStep/2) / m.xStep), collision.Threshold, string(collision.Type)}
							particlePtr.eKinetic -= collision.Threshold
							particlePtr.redirect(cosChi, math.Cos(phi), m)

						case lxgata.IONIZATION:
							collflow <- CollisionEvent{int((particlePtr.x + m.xStep/2) / m.xStep), collision.Threshold, string(collision.Type)}
							ejected := *particlePtr

							collisionEnergy := particlePtr.eKinetic - collision.Threshold

							eScattered := collisionEnergy * rand.Float64()
							particlePtr.eKinetic = eScattered
							cosChiScattered := math.Sqrt(eScattered / collisionEnergy)
							particlePtr.redirect(cosChiScattered, math.Cos(phi), m)

							eEjected := collisionEnergy - eScattered
							ejected.eKinetic = eEjected
							cosChiEjected := math.Sqrt(eEjected / collisionEnergy)
							ejected.redirect(cosChiEjected, math.Cos(phi+math.Pi), m)

							computeflow <- &ejected
							computeWg.Add(1)
						}
					} else {
						// fmt.Printf("null x:%f; cell: %d of %d; mu: %f; eKinRad: %f eKinAx %f\n", particlePtr.x, int(particlePtr.x/m.xStep), m.numCells, particlePtr.mu, particlePtr.eStar, particlePtr.eKinetic*particlePtr.mu)
						currentCell := int(particlePtr.x / m.xStep)
						if currentCell < m.numCells {
							collflow <- CollisionEvent{int(particlePtr.x / m.xStep), 0, "NULL"}
						}
					}
				}
				computeWg.Done()
			}
		}()
	}
	computeWg.Wait()
	close(stateflow)
	close(computeflow)
	close(collflow)
	close(probflow)
	stateWg.Wait()
	print("\r")
}
