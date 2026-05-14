package model

import (
	"fmt"
	"math"
	"math/rand/v2"

	"github.com/wildstyl3r/stmc/internal/config"
	"github.com/wildstyl3r/stmc/internal/utils"
)

type ParticleType int8

const (
	Electron ParticleType = iota
	Ion
)

type Particle struct {
	potential      float64
	eKinetic       float64 // [eV]
	mu             float64
	cosEta, sinEta float64

	prevAxialEnergy       float64
	prevMuSign            float64
	trajectory            TrajectoryConstants
	ejectedFromIonization bool
	producedIonization    bool

	timeToWall float64
	weight     int

	generation int

	ptype ParticleType

	origin int
}

func (m *Model) newParticle(origin int) Particle {
	eKinetic := 4. + rand.Float64()
	var mu float64
	switch m.Parameters.GetEmissionMode() {
	case config.Cosine:
		// mu = math.Cos(math.Asin(math.Sqrt(rand.Float64())))
		mu = math.Sqrt(1 - rand.Float64())
	case config.Forward:
		mu = 1.
	case config.ForwardIsotropic:
		mu = rand.Float64()
	default:
		panic("unexpected config.EmissionMode")
	}
	p := Particle{
		// x:          0,
		potential:  m.VfromL(0),
		eKinetic:   eKinetic,
		mu:         mu,
		prevMuSign: mu,
		origin:     origin,
		weight:     1,
	}
	if m.Parameters.Volumetric {
		// p.y, p.z = utils.UniformOnDisk(m.Parameters.CathodeRadius)

		eta := rand.Float64() * 2. * math.Pi
		p.sinEta, p.cosEta = math.Sin(eta), math.Cos(eta)
	}
	p.recalcParams(m)
	return p
}

func (p *Particle) getAxialEnergy() float64 {
	return p.eKinetic - p.trajectory.radialEnergy
}

func (p *Particle) getPotentialEnergy() float64 {
	return p.trajectory.totEnergy - p.eKinetic
}

func (p *Particle) axialVelocityNow() (v float64) {
	return utils.EV2electronVelocity(p.getAxialEnergy())
}

func (p *Particle) setEnergy(eKinetic float64, s *Model, zeroChangeAcceptable bool) {
	// p.MoveRadial(eKinetic, s, zeroChangeAcceptable)
	// r2 := p.y*p.y + p.z*p.z
	if eKinetic < p.trajectory.radialEnergy {
		fmt.Printf("stopPoint = %f\np: %v; %p\n", p.trajectory.getTurnaroundX(s), p, p)
		panic("eKinetic < p.eStar")
	}
	if !zeroChangeAcceptable && math.Abs(p.eKinetic-eKinetic) < 1e-16 {
		fmt.Printf("need to be at cell: %f coord by V is %f\n", s.LfromV(-p.getPotentialEnergy())/s.XStep, s.LfromV(-p.getPotentialEnergy()))
		panic("no change in energy")
	}

	if eKinetic < 0 {
		panic("eKinetic is < 0")
	}

	p.eKinetic = eKinetic

	if math.IsNaN(eKinetic) {
		panic("eKinetic is NaN")
	}
	p.mu = math.Copysign(math.Sqrt(p.getAxialEnergy()/eKinetic), p.mu)
	p.potential = -p.getPotentialEnergy()
	if s.Parameters.Volumetric {
		p.updateExtraDims(s)
	}
}

func (p *Particle) recalcParams(s *Model) {
	if math.IsNaN(p.mu) {
		// panic("mu is nan")
		fmt.Println("particle's mu is nan")
		p.trajectory.totEnergy = 0
		return
	}
	p.trajectory.radialEnergy = p.eKinetic * (1 - p.mu*p.mu)
	if p.trajectory.radialEnergy < 0 {
		panic("radial energy < 0")
	}
	// r2 := p.y*p.y + p.z*p.z
	p.trajectory.totEnergy = p.eKinetic - p.potential //s.VfromL(p.x)
	if p.trajectory.totEnergy < 0. {
		panic("total energy below 0")
	}
	if p.trajectory.totEnergy > -s.VfromL(0)+5 && !s.Parameters.HasSuperelastics() {
		println("total energy greater than might ever be")
	}

	if s.Parameters.Volumetric {
		p.prevAxialEnergy = p.getAxialEnergy()
		p.prevMuSign = p.mu
	}
}

func (p *Particle) redirect(cosChi, cosPhi float64, m *Model) {
	// p.mu = p.mu*cosChi - math.Sqrt((1.-p.mu*p.mu)*(1.-cosChi*cosChi))*cosPhi // Euler angles formula
	//cos (theta1) = cos(theta)*cos(chi) - sin(theta)*cos(phi)*sin(chi)
	sinChi := math.Sqrt(math.FMA(cosChi, -cosChi, 1.))
	sinTheta := math.Sqrt(math.FMA(p.mu, -p.mu, 1.))
	if math.IsNaN(sinTheta) {
		println("sin(theta) is nan")
	}
	cosPhi_sinChi := cosPhi * sinChi
	oldCosTheta := p.mu

	p.mu = math.FMA(p.mu, cosChi, -sinTheta*cosPhi_sinChi)
	if m.Parameters.Volumetric {
		sinTheta_cosChi_plus_cosTheta_cosPhi_sinPhi := math.FMA(sinTheta, cosChi, oldCosTheta*cosPhi_sinChi)
		sinPhi_sinChi := math.Sqrt(math.FMA(cosPhi, -cosPhi, 1.) * math.FMA(cosChi, -cosChi, 1.))
		nu := math.FMA(p.cosEta, sinTheta_cosChi_plus_cosTheta_cosPhi_sinPhi, p.sinEta*sinPhi_sinChi)
		xi := math.FMA(p.sinEta, sinTheta_cosChi_plus_cosTheta_cosPhi_sinPhi, p.cosEta*sinPhi_sinChi)

		sinThetaPrime := math.Sqrt(math.FMA(p.mu, -p.mu, 1.))
		p.cosEta = nu / sinThetaPrime
		p.sinEta = xi / sinThetaPrime
		e_norm := math.FMA(p.cosEta, p.cosEta, p.sinEta*p.sinEta)
		if math.Abs(e_norm-1.) > 1e-5 {
			p.cosEta /= e_norm
			p.sinEta /= e_norm
		}
	}
	p.recalcParams(m)
}

func (p *Particle) updateExtraDims(m *Model) {

}
