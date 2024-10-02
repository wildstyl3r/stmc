package main

import (
	"fmt"
	"math"
	"math/rand"
	"sync"

	"github.com/wildstyl3r/lxgata"
)

type Model struct {
	crossSections     lxgata.Collisions
	cathodeFallLength float64 // cathode fall length
	Vc                float64 // cathode fall potential
	gapLength         float64 // gap length
	EConst            float64
	density           float64 // N
	nElectrons        int
	threads           int

	Va float64 // additional voltage to avoid numeric negative energy beyond cathode fall

	numCells   int
	numCellsE  int
	numCellsMu int

	xStep  float64
	eStep  float64
	muStep float64

	distribution [][][]int

	ionizationAtCell []int
	elasticAtCell    []int
	excitationAtCell []int
	nullAtCell       []int
	electronsAtCell  []int

	lookUpVelocity  []float64
	lookUpPotential []float64

	energyLossByProcess [][]float64
	probabilities       []map[string][]float64
}

func (m *Model) init() {
	m.Va = -m.EConst*(m.gapLength-m.cathodeFallLength) + 1.

	meanFreePath := 1. / (m.crossSections.SurplusCrossSection() * m.density)
	fmt.Printf("Mean free path: %f\n", meanFreePath)

	m.numCells = 5 * int(m.gapLength/meanFreePath)
	m.numCellsE = int(200./m.eStep) + 1
	m.numCellsMu = int(2./m.muStep) + 1

	m.xStep = m.gapLength / float64(m.numCells)

	m.distribution = make([][][]int, m.numCells)
	for i := range m.distribution {
		m.distribution[i] = make([][]int, m.numCellsE)
		for j := range m.distribution[i] {
			m.distribution[i][j] = make([]int, m.numCellsMu)
		}
	}

	m.ionizationAtCell = make([]int, m.numCells)
	m.nullAtCell = make([]int, m.numCells)
	m.excitationAtCell = make([]int, m.numCells)
	m.elasticAtCell = make([]int, m.numCells)
	m.electronsAtCell = make([]int, m.numCells)
	m.energyLossByProcess = make([][]float64, m.numCells)
	for i := range m.energyLossByProcess {
		m.energyLossByProcess[i] = make([]float64, 3)
	}
	m.probabilities = make([]map[string][]float64, m.numCells)
	for i := range m.probabilities {
		m.probabilities[i] = make(map[string][]float64)
	}

	var energyRoot2Velocity float64 = math.Sqrt(2. / me)
	m.lookUpVelocity = make([]float64, m.numCellsE)
	for i := range m.lookUpVelocity {
		m.lookUpVelocity[i] = math.Sqrt(eV2J(m.eStep*float64(i))) * energyRoot2Velocity
	}
	m.lookUpPotential = make([]float64, m.numCellsE)
	for i := range m.lookUpPotential {
		m.lookUpPotential[i] = m.VfromL(float64(i) * m.xStep)
	}
}

// g
func (s *Model) VfromL(l float64) (V float64) {
	cathodeFallPortion := l / s.cathodeFallLength
	if l < s.cathodeFallLength {
		V = math.Abs(s.Vc) * (cathodeFallPortion*(2.-cathodeFallPortion) - 1.)
	} else {
		V = -s.EConst * (l - s.cathodeFallLength)
	}
	V -= math.Abs(s.Va)
	return V
}

// g^-1
func (s *Model) LfromV(V float64) (l float64) {
	Vtot := math.Abs(s.Va) + math.Abs(s.Vc)
	if V < -math.Abs(s.Va) {
		l = s.cathodeFallLength * (1. - math.Sqrt(1.-(V+Vtot)/math.Abs(s.Vc)))
	} else {
		l = (V+math.Abs(s.Va))/(-s.EConst) + s.cathodeFallLength
	}
	if math.IsNaN(l) {
		panic("x is NaN")
	}
	return l
}

func (s *Model) EFieldFromL(l float64) (E float64) { // V/m
	E = s.EConst
	if l < s.cathodeFallLength {
		E += -2. * (1. - l/s.cathodeFallLength) * math.Abs(s.Vc) / s.cathodeFallLength
	}
	return E
}

type StateIncrement struct {
	xIndex, eIndex, muIndex int
}

type ProbIncrement struct {
	x        int
	collType string
	prob     float64
}

func (s *Model) collisionSelector(p *Particle, eColl, M float64, probChan chan<- ProbIncrement) *lxgata.Collision {
	var totalCrossSectionPrimed float64 = M * math.Abs(s.EFieldFromL(p.x)) / (s.density * math.Sqrt(eColl)) // [m^2]
	var selector float64 = totalCrossSectionPrimed - s.crossSections.TotalCrossSectionAt(eColl)             // the difference is null collision
	var choice float64 = rand.Float64()
	prob := selector / totalCrossSectionPrimed

	currentCell := int(p.x / s.xStep)
	eAtCellBound := s.lookUpPotential[currentCell] + p.totEnergy
	var totalCrossSectionAtCellBound float64 = s.crossSections.TotalCrossSectionAt(eAtCellBound)
	var totalCrossSectionPrimedAtCellBound float64 = M * math.Abs(s.EFieldFromL(p.x)) / (s.density * math.Sqrt(eAtCellBound)) // [m^2]
	probChan <- ProbIncrement{currentCell, "NULL", (totalCrossSectionPrimedAtCellBound - totalCrossSectionAtCellBound) / totalCrossSectionAtCellBound}

	var selection *lxgata.Collision
	selected := false
	if choice < prob {
		selection = nil // null collision
		selected = true
	}
	for i, collision := range s.crossSections {
		selector += collision.CrossSectionAt(eColl)
		prob = selector / totalCrossSectionPrimed
		probChan <- ProbIncrement{currentCell, string(collision.Type), collision.CrossSectionAt(eAtCellBound) / totalCrossSectionAtCellBound}
		if choice < prob && !selected {
			selection = &s.crossSections[i]
			selected = true
		}
	}
	return selection
}

func (s *Model) nextCollision(p *Particle, acc map[StateIncrement]int, probChan chan<- ProbIncrement) *lxgata.Collision { //ch chan<- StateIncrement,
	R := -math.Log(1. - rand.Float64())

	var eColl float64
	if p.mu < 0 {
		stopPoint := s.LfromV(p.eStar - p.totEnergy) // the point where electron speed along x axis becomes 0, hereby mu = 0
		if stopPoint < 0. {
			stopPoint = 0.
		}
		for cellIndex := int(p.x / s.xStep); 0 <= cellIndex && stopPoint < s.xStep*float64(cellIndex); cellIndex-- {
			eLeft := s.lookUpPotential[cellIndex] + p.totEnergy
			if eLeft < 0. {
				panic("eleft < 0")
			}
			M := p.M(cellIndex, s)
			G := 2. * M * (math.Sqrt(p.e-p.eStar) - math.Sqrt(eLeft-p.eStar))
			if G < R { // no collision
				p.setEnergy(eLeft, s)
				eIndex := int((p.e + s.eStep/2.) / s.eStep)
				muIndex := int((p.mu + 1. + s.muStep/2.) / s.muStep)
				if cellIndex < s.numCells && eIndex < s.numCellsE && muIndex < s.numCellsMu {
					acc[StateIncrement{
						xIndex:  cellIndex,
						eIndex:  eIndex,
						muIndex: muIndex,
					}]++
				}

				R -= G
			} else { // collision occured (possibly null)
				eColl = math.Pow(math.Sqrt(p.e-p.eStar)-R/(2.*M), 2.) + p.eStar
				if eColl <= 0 || p.x < 0 {
					p.setEnergy(p.eStar, s)
					break
				}

				p.setEnergy(eColl, s)

				return s.collisionSelector(p, eColl, M, probChan)
			}
		}
	}
	p.mu = math.Copysign(p.mu, +1.)
	// so mu is >= 0
	for cellIndex := int(p.x / s.xStep); cellIndex < s.numCells; cellIndex++ {
		eRight := s.lookUpPotential[cellIndex+1] + p.totEnergy
		M := p.M(cellIndex, s)
		G := 2. * M * (math.Sqrt(eRight-p.eStar) - math.Sqrt(p.e-p.eStar))
		if G < R { // no collision
			p.setEnergy(eRight, s)
			eIndex := int((p.e + s.eStep/2.) / s.eStep)
			muIndex := int((p.mu + 1. + s.muStep/2.) / s.muStep)
			if cellIndex+1 < s.numCells && eIndex < s.numCellsE && muIndex < s.numCellsMu {
				acc[StateIncrement{
					xIndex:  cellIndex + 1,
					eIndex:  eIndex,
					muIndex: muIndex,
				}]++
			}

			R -= G
		} else { // collision occured (possibly null)
			eColl = math.Pow(R/(2.*M)+math.Sqrt(p.e-p.eStar), 2.) + p.eStar

			p.setEnergy(eColl, s)
			return s.collisionSelector(p, eColl, M, probChan)
		}
	}
	return nil
}

type CollisionEvent struct {
	x          int
	energyLoss float64
}

func (s *Model) run() {
	var computeWg, stateWg sync.WaitGroup
	stateflow := make(chan map[StateIncrement]int, 1000)
	ionflow := make(chan CollisionEvent, 100000)
	elasticflow := make(chan CollisionEvent, 100000)
	exitationflow := make(chan CollisionEvent, 100000)
	nullflow := make(chan int, 100000)

	stateWg.Add(1)
	go func() {
		for accum := range stateflow {
			for update, value := range accum {
				s.distribution[update.xIndex][update.eIndex][update.muIndex] += value
				if update.muIndex < s.numCellsMu/2 {
					s.electronsAtCell[update.xIndex]--
				} else {
					s.electronsAtCell[update.xIndex]++
				}
			}
		}
		stateWg.Done()
	}()

	probflow := make(chan ProbIncrement, 100000)
	stateWg.Add(1)
	go func() {
		for probInc := range probflow {
			s.probabilities[probInc.x][probInc.collType] = append(s.probabilities[probInc.x][probInc.collType], probInc.prob)
		}
		stateWg.Done()
	}()
	stateWg.Add(1)
	go func() {
		for ionization := range ionflow {
			s.ionizationAtCell[ionization.x]++
			s.energyLossByProcess[ionization.x][0] += ionization.energyLoss
		}
		stateWg.Done()
	}()
	stateWg.Add(1)
	go func() {
		for elastic := range elasticflow {
			s.elasticAtCell[elastic.x]++
			s.energyLossByProcess[elastic.x][1] += elastic.energyLoss
		}
		stateWg.Done()
	}()
	stateWg.Add(1)
	go func() {
		for excitation := range exitationflow {
			s.excitationAtCell[excitation.x]++
			s.energyLossByProcess[excitation.x][2] += excitation.energyLoss
		}
		stateWg.Done()
	}()
	stateWg.Add(1)
	go func() {
		for nullCell := range nullflow {
			s.nullAtCell[nullCell]++
		}
		stateWg.Done()
	}()

	computeflow := make(chan *Particle, s.nElectrons*100)
	for i := 0; i < s.nElectrons; i++ {
		particle := newParticle()
		particle.recalcParams(s)
		s.distribution[0][int((particle.e+s.eStep/2.)/s.eStep)][int((particle.mu+1.+s.muStep/2.)/s.muStep)]++
		computeflow <- &particle
		computeWg.Add(1)
	}

	status := []string{"//", "==", "\\\\", "||"}
	for w := 0; w < s.threads; w++ {
		go func() {
			counter := 0
			for particlePtr := range computeflow {
				counter++
				print("\r" + status[counter&0b11])
				accumulator := make(map[StateIncrement]int)
				for int(particlePtr.x/s.xStep)+1 < s.numCells {
					if collision := s.nextCollision(particlePtr, accumulator, probflow); collision != nil {
						cosChi := 1. - 2.*rand.Float64()
						cosPhi := math.Cos(2. * math.Pi * rand.Float64())
						switch collision.Type {
						case lxgata.ELASTIC:
							elasticflow <- CollisionEvent{int(particlePtr.x / s.xStep), (1. - (2.*collision.MassRatio)*(1.-cosChi))}
							particlePtr.redirect(cosChi, cosPhi)
							particlePtr.e = particlePtr.e * (1. - (2.*collision.MassRatio)*(1.-cosChi))

						case lxgata.EFFECTIVE:
							elasticflow <- CollisionEvent{int(particlePtr.x / s.xStep), (1. - (2.*collision.MassRatio)*(1.-cosChi))}
							particlePtr.redirect(cosChi, cosPhi)
							particlePtr.e = particlePtr.e * (1. - (2.*collision.MassRatio)*(1.-cosChi))

						case lxgata.EXCITATION:
							exitationflow <- CollisionEvent{int(particlePtr.x / s.xStep), collision.Threshold}
							particlePtr.redirect(cosChi, cosPhi)
							particlePtr.e -= collision.Threshold

						case lxgata.IONIZATION:
							ionflow <- CollisionEvent{int(particlePtr.x / s.xStep), collision.Threshold}
							collisionEnergy := particlePtr.e
							collisionEnergy -= collision.Threshold
							e1 := collisionEnergy * rand.Float64()
							e2 := collisionEnergy - e1
							cosChi1 := math.Sqrt(e1 / collisionEnergy)
							cosChi2 := math.Sqrt(e2 / collisionEnergy)

							ejected := *particlePtr
							ejected.e = e2
							ejected.redirect(cosChi2, cosPhi)
							ejected.recalcParams(s)
							computeflow <- &ejected
							computeWg.Add(1)

							particlePtr.e = e1
							particlePtr.redirect(cosChi1, cosPhi)
						default:
						}
						particlePtr.recalcParams(s)
					} else {
						currentCell := int(particlePtr.x / s.xStep)
						if currentCell < s.numCells {
							nullflow <- currentCell
						}
					}
				}
				stateflow <- accumulator
				computeWg.Done()
			}
		}()
	}
	computeWg.Wait()
	close(stateflow)
	close(computeflow)
	close(ionflow)
	close(nullflow)
	close(exitationflow)
	close(elasticflow)
	close(probflow)
	stateWg.Wait()
	print("\r")
}
