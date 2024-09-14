package main

import (
	"fmt"
	"math"
	"math/rand"
)

type Particle struct {
	x  float64 //  [m]
	e  float64 // [eV]
	mu float64

	// params
	eStar     float64 //[eV]
	totEnergy float64 //aka Vcap
}

func newParticle() Particle {
	initEnergy := 4. + rand.Float64()
	mu := 1. - 2.*rand.Float64()
	eStar := initEnergy * (1 - mu*mu)
	return Particle{
		0,
		initEnergy,
		mu,

		eStar,
		0,
	}
}

func (p *Particle) setEnergy(e float64, s *Model) {
	if e <= p.eStar {
		fmt.Printf("stopPoint = %f\np: %v", s.LfromV(p.eStar-p.totEnergy), p)
		p.e = p.eStar
		p.mu = 0
	} else {
		p.e = e
		p.mu = math.Copysign(math.Sqrt((e-p.eStar)/e), p.mu)
	}
	p.x = s.LfromV(p.e - p.totEnergy)
}

func (p *Particle) recalcParams(s *Model) {
	if math.IsNaN(p.mu) {
		panic("mu is nan")
	}
	p.eStar = p.e * (1 - p.mu*p.mu)
	p.totEnergy = p.e - s.VfromL(p.x)
	if p.totEnergy < 0. {
		panic("tot energy below 0")
	}
}

func (p *Particle) M(i int, s *Model) float64 {
	left, right := s.lookUpPotential[i]+p.totEnergy, s.lookUpPotential[i+1]+p.totEnergy
	if left < 0. {
		left = 0.
	}
	if right <= 0. {
		panic("stf")
	}
	return ternarySearchMaxF(func(e float64) float64 {
		return s.density * s.crossSections.TotalCrossSectionAt(e) * math.Sqrt(e) / math.Abs(s.EFieldFromL(s.LfromV(e-p.totEnergy)))
	}, s.lookUpPotential[i]+p.totEnergy, s.lookUpPotential[i+1]+p.totEnergy, 0.00001)
}

func (p *Particle) redirect(cosChi, cosPhi float64) {
	p.mu = p.mu*cosChi + math.Sqrt((1.-p.mu*p.mu)*(1.-cosChi*cosChi))*cosPhi
}
