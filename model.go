package main

import (
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

	Va float64 // additional voltage to avoid numeric negative energy beyond cathode fall

	numCells   int
	numCellsE  int
	numCellsMu int

	xStep  float64
	eStep  float64
	muStep float64

	flux          []int // ratio to 1 electron emitted from cathode
	driftVelocity []float64
	TownsendAlpha []float64
	// psi             [][][]int
	// psiF            [][][]float64
	psiFIncrement   float64
	distribution    [][][]float64
	electronDensity []float64
	electronEnergy  []float64

	ionizationAtCell []int
	sourceTerm       []float64
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

	// m.psi = make([][][]int, m.numCells)
	m.distribution = make([][][]float64, m.numCells)
	// m.psiF = make([][][]float64, m.numCells)
	m.psiFIncrement = m.cathodeFlux / (float64(m.nElectrons) * m.eStep * m.muStep)
	for i := range m.distribution {
		// m.psi[i] = make([][]int, m.numCellsE)
		m.distribution[i] = make([][]float64, m.numCellsE)
		// m.psiF[i] = make([][]float64, m.numCellsE)
		for j := range m.distribution[i] {
			// m.psi[i][j] = make([]int, m.numCellsMu)
			m.distribution[i][j] = make([]float64, m.numCellsMu)
			// m.psiF[i][j] = make([]float64, m.numCellsMu)
		}
	}
	m.electronDensity = make([]float64, m.numCells)
	m.electronEnergy = make([]float64, m.numCells)

	m.ionizationAtCell = make([]int, m.numCells)
	m.sourceTerm = make([]float64, m.numCells)

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
	if 0 <= xIndex && xIndex < s.numCells {
		if eIndex < s.numCellsE && muIndex < s.numCellsMu && eIndex >= 0 && muIndex >= 0 {
			vxAbs := math.Sqrt(eV2J(p.e)*2./me) * math.Abs(p.mu)
			s.distribution[xIndex][eIndex][muIndex] += s.psiFIncrement / vxAbs
		}
		if p.mu < 0 {
			s.flux[xIndex]--
		} else {
			s.flux[xIndex]++
		}
	}

}

func (s *Model) nextCollision(p *Particle) *lxgata.Collision {
	R := -math.Log(1. - rand.Float64())

	var eColl float64
	if p.mu < 0 {
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

				currentCellIndex := int(p.x / s.xStep)

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
							s.ionizationAtCell[currentCellIndex]++
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
		eRight := s.VfromL(float64(cellIndex+1)*s.xStep) + p.totEnergy
		M := p.M(cellIndex, s)
		G := 2. * M * (math.Sqrt(eRight-p.eStar) - math.Sqrt(p.e-p.eStar))
		if G < R { // no collision
			p.setEnergy(eRight, s)
			s.incrementPsi(p, cellIndex+1)

			R -= G
		} else { // collision occured (possibly null)
			eColl = math.Pow(R/(2.*M)+math.Sqrt(p.e-p.eStar), 2.) + p.eStar

			p.setEnergy(eColl, s)

			currentCellIndex := int(p.x / s.xStep)

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
						s.ionizationAtCell[currentCellIndex]++
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

	particlePtrs := make(map[*Particle]struct{})
	for i := 0; i < s.nElectrons; i++ {
		particle := newParticle()
		particle.recalcParams(s)
		s.incrementPsi(&particle, 0)
		particlePtrs[&particle] = struct{}{}
	}

	for len(particlePtrs) > 0 {
		for particlePtr := range particlePtrs {
			if int(particlePtr.x/s.xStep)+1 >= s.numCells {
				delete(particlePtrs, particlePtr)
				continue
			}

			if collision := s.nextCollision(particlePtr); collision != nil {
				cosChi := 1. - 2.*rand.Float64()
				cosPhi := math.Cos(2. * math.Pi * rand.Float64())
				switch collision.Type {
				case lxgata.ELASTIC:
					particlePtr.redirect(cosChi, cosPhi)
					particlePtr.e = particlePtr.e * (1. - (2.*collision.MassRatio)*(1.-cosChi))

				case lxgata.EFFECTIVE:
					particlePtr.redirect(cosChi, cosPhi)
					particlePtr.e = particlePtr.e * (1. - (2.*collision.MassRatio)*(1.-cosChi))

				case lxgata.EXCITATION:
					particlePtr.redirect(cosChi, cosPhi)
					particlePtr.e -= collision.Threshold

				case lxgata.IONIZATION:
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

					particlePtr.e = e1
					particlePtr.redirect(cosChi1, cosPhi)
				default:
				}
				particlePtr.recalcParams(s)
			}
		}
	}

	var energyRoot2Velocity float64 = math.Sqrt(2. / me)

	for xIndex := 0; xIndex < s.numCells; xIndex++ {
		for eIndex := 1; eIndex < s.numCellsE; eIndex++ {
			currentEnergy := s.eStep * float64(eIndex)
			for muIndex := 0; muIndex < s.numCellsMu; muIndex++ {
				currentMu := s.muStep*float64(muIndex) - 1.

				s.electronDensity[xIndex] += s.distribution[xIndex][eIndex][muIndex]

				if eIndex > 0 {
					s.electronEnergy[xIndex] += s.distribution[xIndex][eIndex][muIndex] * currentEnergy

					s.driftVelocity[xIndex] += s.distribution[xIndex][eIndex][muIndex] * math.Sqrt(eV2J(currentEnergy)) * currentMu

					s.TownsendAlpha[xIndex] += s.crossSections.TotalCrossSectionOfKindAt(lxgata.IONIZATION, currentEnergy) * math.Sqrt(eV2J(currentEnergy)) * s.distribution[xIndex][eIndex][muIndex] * s.eStep
				}
			}
		}
		s.electronEnergy[xIndex] /= s.electronDensity[xIndex]

		s.driftVelocity[xIndex] *= energyRoot2Velocity
		s.driftVelocity[xIndex] /= s.electronDensity[xIndex]

		s.electronDensity[xIndex] *= s.eStep

		s.TownsendAlpha[xIndex] *= energyRoot2Velocity

		s.TownsendAlpha[xIndex] /= s.driftVelocity[xIndex]
		//s.sourceTerm[xIndex] = s.TownsendAlpha[xIndex] * psiX / float64(s.nElectrons)
		s.sourceTerm[xIndex] = s.TownsendAlpha[xIndex] * float64(s.flux[xIndex]) / float64(s.nElectrons)
		// if xIndex+1 < s.numCells {
		// 	s.TownsendAlpha[xIndex] = (float64(s.flux[xIndex+1]) - float64(s.nElectrons)) / (float64(s.nElectrons) * s.xStep * float64(xIndex))
		// 	s.sourceTerm[xIndex] = (float64(s.flux[xIndex+1]) - float64(s.nElectrons)) * float64(s.nElectrons) / (float64(s.nElectrons) * float64(s.nElectrons) * s.xStep * float64(xIndex+1)) //s.TownsendAlpha[xIndex] * s.flux[xIndex] / (s.cathodeFlux)
		// }

		//s.flux[xIndex] = s.driftVelocity[xIndex] * s.electronDensity[xIndex]

	}
}
