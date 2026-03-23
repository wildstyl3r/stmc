package model

import (
	"fmt"
	"math"
	"math/rand"

	"github.com/wildstyl3r/stmc/internal/config"
	"github.com/wildstyl3r/stmc/internal/utils"
)

type Particle struct {
	x, y, z        float64 //  [m]
	eKinetic       float64 // [eV]
	mu             float64
	cosEta, sinEta float64

	prevAxialEnergy float64
	prevX           float64
	prevMuSign      float64
	// params
	trajectory            TrajectoryConstants
	ejectedFromIonization bool
	producedIonization    bool

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
		x:          0,
		eKinetic:   eKinetic,
		mu:         mu,
		prevMuSign: mu,
		origin:     origin,
	}
	if m.Parameters.Volumetric {
		p.y, p.z = utils.UniformOnDisk(m.Parameters.CathodeRadius)

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

func (p *Particle) setEnergy(eKinetic float64, s *Model, zeroChangeAcceptable bool, setX bool) {
	// p.MoveRadial(eKinetic, s, zeroChangeAcceptable)
	// r2 := p.y*p.y + p.z*p.z
	if eKinetic < p.trajectory.radialEnergy {
		fmt.Printf("stopPoint = %f\np: %v; %p\n", p.trajectory.getTurnaroundX(s), p, p)
		panic("eKinetic < p.eStar")
	}
	if math.Abs(p.eKinetic-eKinetic) < 1e-16 && !zeroChangeAcceptable {
		fmt.Printf("need to be at cell: %f coord by V is %f, coord real is %f\n", s.LfromV(-p.getPotentialEnergy())/s.XStep, s.LfromV(-p.getPotentialEnergy()), p.x)
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
	if setX {
		p.x = s.LfromV(-p.getPotentialEnergy())
	}
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
	// r2 := p.y*p.y + p.z*p.z
	p.trajectory.totEnergy = p.eKinetic + -s.VfromL(p.x)
	if p.trajectory.totEnergy < 0. {
		panic("total energy below 0")
	}
	if p.trajectory.totEnergy > -s.VfromL(0)+5 && !s.Parameters.HasSuperelastics() {
		println("total energy greater than might ever be")
	}

	if s.Parameters.Volumetric {
		p.prevAxialEnergy = p.getAxialEnergy()
		p.prevX = p.x
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
	/// updates particle's y and z, and sets new prev_bt value
	var timeIntervalsToSum []float64
	if math.Signbit(p.prevMuSign) != math.Signbit(p.mu) { //reversal occured
		// xRev := m.LfromV(p.getStopPotential())

		// timeIntervalsToSum = p.getTimeIntervalsBetweenPositionsNoReversal(p.prevX, xRev, p.prevAxialEnergy, 0, m)
		// timeIntervalsToSum = append(timeIntervalsToSum, p.getTimeIntervalsBetweenPositionsNoReversal(xRev, p.x, 0, p.getAxialEnergy(), m)...)

	} else {
		// timeIntervalsToSum = p.getTimeIntervalsBetweenPositionsNoReversal(p.prevX, p.x, p.prevAxialEnergy, p.getAxialEnergy(), m)
	}

	timeBetweenCollisions := utils.SumFloat64Slice(timeIntervalsToSum, true)

	p.y = p.sinEta * utils.EV2electronVelocity(p.trajectory.radialEnergy) * timeBetweenCollisions
	p.z = p.cosEta * utils.EV2electronVelocity(p.trajectory.radialEnergy) * timeBetweenCollisions
}
