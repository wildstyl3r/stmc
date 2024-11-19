package main

import (
	"fmt"
	"math"
	"math/rand"
)

type Particle struct {
	x        float64 //  [m]
	eKinetic float64 // [eV]
	mu       float64

	// params
	eStar     float64 //[eV]
	totEnergy float64 //aka Vcap
}

func (m *Model) newParticle() Particle {
	p := Particle{
		x:        0,
		eKinetic: 4. + rand.Float64(),
		mu:       rand.Float64(),
	}
	p.recalcParams(m)
	return p
}

func (p *Particle) setEnergy(eKinetic float64, s *Model, zeroChangeAcceptable bool) {
	if eKinetic < p.eStar {
		stopPotential := -(p.totEnergy - p.eStar)
		fmt.Printf("stopPoint = %f\np: %v; %p\n", s.LfromV(stopPotential), p, p)
		panic("eKinetic < p.eStar")
	}
	if math.Abs(p.eKinetic-eKinetic) < 1e-16 && !zeroChangeAcceptable {
		fmt.Printf("need to be at cell: %f coord by V is %f, coord real is %f\n", s.LfromV(-(p.totEnergy-p.eKinetic))/s.xStep, s.LfromV(-(p.totEnergy - p.eKinetic)), p.x)
		panic("no change in energy")
	}
	p.eKinetic = eKinetic
	p.mu = math.Copysign(math.Sqrt((eKinetic-p.eStar)/eKinetic), p.mu)
	potential := -(p.totEnergy - p.eKinetic)
	p.x = s.LfromV(potential)
}

func (p *Particle) recalcParams(s *Model) {
	if math.IsNaN(p.mu) {
		panic("mu is nan")
	}
	p.eStar = p.eKinetic * (1 - p.mu*p.mu)
	p.totEnergy = p.eKinetic + -s.VfromL(p.x)
	if p.totEnergy < 0. {
		panic("tot energy below 0")
	}
}

func (p *Particle) M(i int, s *Model) float64 {
	eLeft, eRight := p.totEnergy+s.lookUpPotential[i], p.totEnergy+s.lookUpPotential[i+1]
	if eLeft < 0. {
		eLeft = 0.
	}
	if eRight <= 0. {
		panic("stf")
	}
	return ternarySearchMaxF(func(eKin float64) float64 {
		potential := -(p.totEnergy - eKin)
		return s.gasDensity * s.crossSections.TotalCrossSectionAt(eKin) * math.Sqrt(eKin) / math.Abs(s.EFieldFromL(s.LfromV(potential)))
	}, eLeft, eRight, 0.00001)
}

func (p *Particle) redirect(cosChi, cosPhi float64, m *Model) {
	p.mu = p.mu*cosChi - math.Sqrt((1.-p.mu*p.mu)*(1.-cosChi*cosChi))*cosPhi // Euler angles formula
	p.recalcParams(m)
}
