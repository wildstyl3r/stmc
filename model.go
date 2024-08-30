package main

import (
	"fmt"
	"math"
	"math/rand"

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
	cathodeFlux       float64

	Va float64 // additional voltage to avoid numeric negative energy problems beyond cathode fall

	numCells   int
	numCellsE  int
	numCellsMu int

	xStep  float64
	eStep  float64
	muStep float64

	flux            []int // ratio to 1 electron emitted from cathode
	driftVelocity   []float64
	TownsendAlpha   []float64
	psi             [][][]int
	distribution    [][][]float64
	electronDensity []float64
	electronEnergy  []float64

	ionizationAtCell               []int
	sourceTerm                     []float64
	accumulatedAvgIonizationAtCell []float64
}

func (m *Model) init() {
	m.Va = -m.EConst*(m.gapLength-m.cathodeFallLength) + 1.

	meanFreePath := 1. / (m.crossSections.SurplusCrossSection() * m.density)

	m.numCells = 3 * int(m.gapLength/meanFreePath)
	m.numCellsE = int(200./m.eStep) + 1
	m.numCellsMu = int(2./m.muStep) + 1

	m.xStep = m.gapLength / float64(m.numCells)

	m.flux = make([]int, m.numCells)
	m.driftVelocity = make([]float64, m.numCells)
	m.TownsendAlpha = make([]float64, m.numCells)

	m.psi = make([][][]int, m.numCells)
	m.distribution = make([][][]float64, m.numCells)
	for i := range m.psi {
		m.psi[i] = make([][]int, m.numCellsE)
		m.distribution[i] = make([][]float64, m.numCellsE)
		for j := range m.psi[i] {
			m.psi[i][j] = make([]int, m.numCellsMu)
			m.distribution[i][j] = make([]float64, m.numCellsMu)
		}
	}
	m.electronDensity = make([]float64, m.numCells)
	m.electronEnergy = make([]float64, m.numCells)

	m.ionizationAtCell = make([]int, m.numCells)
	m.sourceTerm = make([]float64, m.numCells)
	m.accumulatedAvgIonizationAtCell = make([]float64, m.numCells)

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

func (s *Model) Q(e float64) float64 {
	return s.density * s.crossSections.TotalCrossSectionAt(e)
}

func (s *Model) incrementPsi(p *Particle, xIndex int) {
	eIndex := int((p.e + s.eStep/2.) / s.eStep)
	muIndex := int((p.mu + 1. + s.muStep/2.) / s.muStep)
	if xIndex < s.numCells && eIndex < s.numCellsE && muIndex < s.numCellsMu && xIndex >= 0 && eIndex >= 0 && muIndex >= 0 {
		if p.mu < 0 {
			s.psi[xIndex][eIndex][muIndex]--
			s.flux[xIndex]--
		} else {
			s.psi[xIndex][eIndex][muIndex]++
			s.flux[xIndex]++
		}
	}

}

func (s *Model) nextCollision(p *Particle) *lxgata.Collision {
	R := -math.Log(1. - rand.Float64())

	var eColl float64
	if p.mu < 0 {
		//fmt.Println("-mu")
		stopPoint := s.LfromV(p.eStar - p.totEnergy) // the point where electron speed along x axis becomes 0, hereby mu = 0
		if stopPoint < 0. {
			stopPoint = 0.
		}
		for cellIndex := int(p.x / s.xStep); 0 <= cellIndex && stopPoint < s.xStep*float64(cellIndex); cellIndex-- {
			eLeft := s.VfromL(float64(cellIndex)*s.xStep) + p.totEnergy
			if eLeft < 0. {
				panic("eleft < 0")
			}
			M := p.M(cellIndex, s)
			G := 2. * M * (math.Sqrt(p.e-p.eStar) - math.Sqrt(eLeft-p.eStar))
			if G < R { // no collision
				p.setEnergy(eLeft, s)
				s.incrementPsi(p, cellIndex)

				R -= G
			} else { // collision occured (possibly null)
				eColl = math.Pow(math.Sqrt(p.e-p.eStar)-R/(2.*M), 2.) + p.eStar
				if eColl <= 0 || p.x < 0 {
					p.setEnergy(p.eStar, s)
					break
				}

				p.setEnergy(eColl, s)
				//fmt.Printf("eColl: %f; eTot: %f\n", eColl, p.totEnergy)

				var totalCrossSectionPrimed float64 = M * math.Abs(s.EFieldFromL(p.x)) / (s.density * math.Sqrt(eColl))
				var selector float64 = totalCrossSectionPrimed - s.crossSections.TotalCrossSectionAt(eColl) // the difference is null collision
				var choice float64 = rand.Float64()
				if choice < selector/totalCrossSectionPrimed {
					return nil // null collision
				}
				for i, collision := range s.crossSections {
					selector += collision.CrossSectionAt(eColl)
					if choice < selector/totalCrossSectionPrimed {
						if s.crossSections[i].Type == lxgata.IONIZATION {
							s.ionizationAtCell[int(p.x/s.xStep)]++
						}
						return &s.crossSections[i]
					}
				}
				break
			}
		}
	}
	p.mu = math.Copysign(p.mu, +1.)
	// so mu is >= 0
	for cellIndex := int(p.x / s.xStep); cellIndex < s.numCells; cellIndex++ {
		//fmt.Printf("p.x = %f; cellIndex = %d\n", p.x, int(p.x/s.xStep))
		eRight := s.VfromL(float64(cellIndex+1)*s.xStep) + p.totEnergy
		M := p.M(cellIndex, s)
		G := 2. * M * (math.Sqrt(eRight-p.eStar) - math.Sqrt(p.e-p.eStar))
		if G < R { // no collision
			if p.x > float64(cellIndex+1)*s.xStep {
				fmt.Println("wtf")
			}
			p.setEnergy(eRight, s)
			s.incrementPsi(p, cellIndex+1)

			R -= G
		} else { // collision occured (possibly null)
			//fmt.Printf("eStar: %f; R: %f, M: %f\n", p.eStar, R, M)
			eColl = math.Pow(R/(2.*M)+math.Sqrt(p.e-p.eStar), 2.) + p.eStar
			if p.x > s.LfromV(eColl-p.totEnergy) {
				fmt.Println("wtf +mu coll")
			}

			p.setEnergy(eColl, s)

			var totalCrossSectionPrimed float64 = M * math.Abs(s.EFieldFromL(p.x)) / (s.density * math.Sqrt(eColl))
			var selector float64 = totalCrossSectionPrimed - s.crossSections.TotalCrossSectionAt(eColl) // the difference is null collision
			var choice float64 = rand.Float64()
			if choice < selector/totalCrossSectionPrimed {
				return nil // null collision
			}
			for i, collision := range s.crossSections {
				selector += collision.CrossSectionAt(eColl)
				if choice < selector/totalCrossSectionPrimed {
					if s.crossSections[i].Type == lxgata.IONIZATION {
						s.ionizationAtCell[int(p.x/s.xStep)]++
					}
					return &s.crossSections[i]
				}
			}
			break
		}
	}
	return nil
}

func (s *Model) run() {
	var particles []Particle
	for i := 0; i < s.nElectrons; i++ {
		particles = append(particles, newParticle())
		particles[len(particles)-1].recalcParams(s)
		s.incrementPsi(&particles[len(particles)-1], 0)
	}

	activeParticles := true
	for activeParticles {
		if math.IsNaN(particles[0].x) {
			panic("nan")
		}
		activeParticles = false
		for i := len(particles) - 1; i >= 0; i-- {
			if int(particles[i].x/s.xStep)+1 < s.numCells {
				//fmt.Println("any activity ", i, " ", particles[i].x, " ", particles[i].x/s.xStep, " ", s.numCells, particles[i].mu, particles[i].e)
				activeParticles = true
			} else {
				continue
			}

			if collision := s.nextCollision(&particles[i]); collision != nil {
				cosChi := 1. - 2.*rand.Float64()
				cosPhi := math.Cos(2. * math.Pi * rand.Float64())
				switch collision.Type {
				case lxgata.ELASTIC:
					particles[i].redirect(cosChi, cosPhi)
					particles[i].e = particles[i].e * (1. - (2.*collision.MassRatio)*(1.-cosChi))

				case lxgata.EXCITATION:
					particles[i].redirect(cosChi, cosPhi)
					particles[i].e -= collision.Threshold

				case lxgata.IONIZATION:
					collisionEnergy := particles[i].e
					collisionEnergy -= collision.Threshold
					e1 := collisionEnergy * rand.Float64()
					e2 := collisionEnergy - e1
					cosChi1 := math.Sqrt(e1 / collisionEnergy)
					cosChi2 := math.Sqrt(e2 / collisionEnergy)

					particles = append(particles, particles[i])
					//fmt.Println("ionization occured; active particles: {}", len(particles))

					particles[len(particles)-1].e = e2
					particles[len(particles)-1].redirect(cosChi2, cosPhi)
					particles[len(particles)-1].recalcParams(s)

					particles[i].e = e1
					particles[i].redirect(cosChi1, cosPhi)
				default:
				}
				particles[i].recalcParams(s)
			}
		}
	}

	const me float64 = 9.1093837139e-31 // [kg]

	for xIndex := 0; xIndex < s.numCells; xIndex++ {

		for eIndex := 1; eIndex < s.numCellsE; eIndex++ {
			//currentE := s.eStep * float64(eIndex)
			v := math.Sqrt(2. * eV2J(s.eStep*float64(eIndex)) / me) // [m/s]
			for muIndex := 0; muIndex < s.numCellsMu; muIndex++ {
				//counter++
				currentMu := s.muStep*float64(muIndex) - 1.
				//preVxSum += math.Abs(currentMu) * math.Sqrt(currentE) * float64(s.psi[xIndex][eIndex][muIndex])

				if math.Abs(currentMu) > 0.00001 {
					s.distribution[xIndex][eIndex][muIndex] = float64(s.psi[xIndex][eIndex][muIndex]) * s.cathodeFlux / (v * math.Abs(currentMu) * s.eStep * s.muStep * float64(s.nElectrons))
				}
				if math.IsNaN(s.distribution[xIndex][eIndex][muIndex]) {
					fmt.Printf("distr nan %d %d %d\n", xIndex, eIndex, muIndex)
				}
			}
		}
	}

	for xIndex := 0; xIndex < s.numCells; xIndex++ {
		NxI, NxINext := 0., 0.
		for eIndex := 1; eIndex < s.numCellsE; eIndex++ {
			currentEnergy := s.eStep * float64(eIndex)
			for muIndex := 0; muIndex < s.numCellsMu; muIndex++ {
				currentMu := s.muStep*float64(muIndex) - 1.

				s.electronDensity[xIndex] += s.distribution[xIndex][eIndex][muIndex] * s.eStep

				if eIndex > 0 {
					s.electronEnergy[xIndex] += s.distribution[xIndex][eIndex][muIndex] * s.eStep * currentEnergy

					s.driftVelocity[xIndex] += s.distribution[xIndex][eIndex][muIndex] * s.eStep * math.Sqrt(eV2J(currentEnergy)) * currentMu

					//s.TownsendAlpha[xIndex] += s.crossSections.TotalCrossSectionOfKindAt(lxgata.IONIZATION, currentEnergy) * math.Sqrt(eV2J(currentEnergy)) * s.distribution[xIndex][eIndex][muIndex] * s.eStep
				}

				NxI += float64(s.psi[xIndex][eIndex][muIndex])
				if xIndex+1 < s.numCells {
					NxINext += float64(s.psi[xIndex+1][eIndex][muIndex])
				}
			}
		}
		s.electronEnergy[xIndex] /= s.electronDensity[xIndex]

		s.driftVelocity[xIndex] *= math.Sqrt(2. / me)
		s.driftVelocity[xIndex] /= s.electronDensity[xIndex]

		//s.TownsendAlpha[xIndex] *= math.Sqrt(2. / me)
		//s.TownsendAlpha[xIndex] /= s.driftVelocity[xIndex]
		if xIndex+1 < s.numCells {
			s.TownsendAlpha[xIndex] = (float64(s.flux[xIndex+1]) - float64(s.nElectrons)) / (float64(s.nElectrons) * s.xStep * float64(xIndex))
			s.sourceTerm[xIndex] = (float64(s.flux[xIndex+1]) - float64(s.nElectrons)) * float64(s.nElectrons) / (float64(s.nElectrons) * float64(s.nElectrons) * s.xStep * float64(xIndex+1)) //s.TownsendAlpha[xIndex] * s.flux[xIndex] / (s.cathodeFlux)
		}

		//s.flux[xIndex] = s.driftVelocity[xIndex] * s.electronDensity[xIndex]

	}
}
