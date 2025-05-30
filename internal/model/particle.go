package model

import (
	"fmt"
	"math"
	"math/rand"
)

type Particle struct {
	x, y, z        float64 //  [m]
	eKinetic       float64 // [eV]
	mu             float64
	cosEta, sinEta float64

	c2, c1     float64
	prevBt     float64
	prevMuSign float64
	// params
	eStar             float64 //[eV]
	totEnergy         float64 //aka Vcap
	_debug_IonEjected bool

	origin int
}

func (m *Model) newParticle(origin int) Particle {
	y, z := UniformOnDisk(m.parameters.CathodeRadius)
	eta := rand.Float64() * 2. * math.Pi
	eKinetic := 4. + rand.Float64()
	mu := rand.Float64()
	p := Particle{
		x:          0,
		y:          y,
		z:          z,
		eKinetic:   eKinetic,
		mu:         mu,
		sinEta:     math.Sin(eta),
		cosEta:     math.Cos(eta),
		c1:         -m.parameters.CathodeFallLength,
		c2:         m.parameters.CathodeFallLength * math.Sqrt(eKinetic*mu*mu/m.parameters.CathodeFallPotential),
		prevMuSign: mu,
		origin:     origin,
	}
	p.recalcParams(m)
	return p
}

func (p *Particle) setEnergy(eKinetic float64, s *Model, zeroChangeAcceptable bool, setX bool) {
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
	if setX {
		potential := -(p.totEnergy - p.eKinetic)
		p.x = s.LfromV(potential)
	}
}

func (p *Particle) recalcParams(s *Model) {
	if math.IsNaN(p.mu) {
		// panic("mu is nan")
		fmt.Println("particle's mu is nan")
		p.totEnergy = 0
		return
	}
	p.eStar = p.eKinetic * (1 - p.mu*p.mu)
	p.totEnergy = p.eKinetic + -s.VfromL(p.x)
	if p.totEnergy < 0. {
		panic("tot energy below 0")
	}
	if p.totEnergy > s.Va+s.Vc+5 {
		println("wtf, tot energy greater than might ever be")
	}

	if s.parameters.Volumetric {
		cosBt, sinBt := math.Cos(p.prevBt), math.Sin(p.prevBt)
		b := math.Sqrt(2.*electronCharge/electornMass*s.parameters.CathodeFallPotential) / s.parameters.CathodeFallLength
		vx := math.Copysign(eV2electronVelocity(p.eKinetic-p.eStar), p.mu)
		currentXminusD := p.x - s.parameters.CathodeFallLength
		lower := math.FMA(b*cosBt, cosBt, sinBt*vx)
		p.c1 = math.FMA(currentXminusD*b, cosBt, -sinBt*vx) / lower
		p.c2 = math.FMA(b*sinBt, currentXminusD, vx*cosBt) / lower
		if math.IsNaN(p.c1) || math.IsNaN(p.c2) {
			panic("c1 or c2 is NaN")
		}
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
		return -s.parameters.gasDensity * s.parameters._crossSections.TotalCrossSectionAt(eKin) * math.Sqrt(eKin) / s.EFieldFromL(s.LfromV(potential))
	}, eLeft, eRight, 0.00001)
}

func (p *Particle) redirect(cosChi, cosPhi float64, m *Model) {
	//cos (theta1) = cos(theta)*cos(chi) - sin(theta)*cos(phi)*sin(chi)
	sinChi := math.Sqrt(math.FMA(cosChi, -cosChi, 1.))
	sinTheta := math.Sqrt(math.FMA(p.mu, -p.mu, 1.))
	if math.IsNaN(sinTheta) {
		println("sin(theta) is nan")
	}
	cosPhi_sinChi := cosPhi * sinChi
	oldCosTheta := p.mu

	p.mu = math.FMA(p.mu, cosChi, -sinTheta*cosPhi_sinChi)
	if m.parameters.Volumetric {
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
		p.prevMuSign = p.mu
	}
	p.recalcParams(m)
}

func (p *Particle) xAnalytic(bt float64, m *Model) float64 {
	return math.FMA(
		p.c2, math.Sin(bt),
		math.FMA(
			p.c1, math.Cos(bt),
			m.inverseCathodeFallLength,
		),
	)
}

func (p *Particle) updateExtraDims(m *Model) {
	/// updates particle's y and z, and sets new prev_bt value
	b := math.Sqrt(2.*electronCharge/electornMass*m.parameters.CathodeFallPotential) / m.parameters.CathodeFallLength
	var bt float64
	btCathodeFallLength := math.Atan(-p.c1 / p.c2)
	if math.IsNaN(btCathodeFallLength) {
		panic("btCathodeFallLength is NaN")
	}
	btReverse := math.Atan(p.c2 / p.c1)
	if math.IsNaN(btReverse) {
		panic("btReverse is NaN")
	}
	for btReverse < p.prevBt {
		btReverse += math.Pi
	}
	//check negative mu case
	if p.prevMuSign > 0 && p.mu > 0 {
		bt = ternarySearchMax(func(bt float64) float64 {
			bt_x := p.xAnalytic(bt, m)
			return -(p.x - bt_x) * (p.x - bt_x)
		}, p.prevBt, btCathodeFallLength, 1e-6)
		if math.IsNaN(bt) {
			panic("bt is NaN in case pMu > 0 && mu > 0")
		}
	} else if p.prevMuSign <= 0 && p.mu > 0 {
		btBeforeReverse := btReverse - p.prevBt
		btAfterReverse := ternarySearchMax(func(bt float64) float64 {
			bt_x := p.xAnalytic(bt, m)
			return -(p.x - bt_x) * (p.x - bt_x)
		}, btReverse, btCathodeFallLength, 1e-6)
		bt = btBeforeReverse + btAfterReverse
		if math.IsNaN(bt) {
			panic("bt is NaN in case pMu < 0 && mu > 0")
		}
	} else if p.prevMuSign <= 0 && p.mu <= 0 {
		bt = ternarySearchMax(func(bt float64) float64 {
			bt_x := p.xAnalytic(bt, m)
			return -(p.x - bt_x) * (p.x - bt_x)
		}, p.prevBt, btReverse, 1e-6)
		if math.IsNaN(bt) {
			panic("bt is NaN in case pMu < 0 && mu < 0")
		}
	} else {
		fmt.Printf("particle: %v", *p)
		panic("should-be-impossible condition: prev mu > 0, current mu < 0")
	}
	t := bt / b
	p.y = p.sinEta * eV2electronVelocity(p.eStar) * t
	p.z = p.cosEta * eV2electronVelocity(p.eStar) * t
	p.prevBt = bt
}
