package model

import (
	"fmt"
	"math"
	"math/rand"

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
	// eStar                 float64 //[eV]
	// totEnergy             float64 //aka Vcap
	ejectedFromIonization bool
	trajectory            TrajectoryConstants

	origin int
}

func (m *Model) newParticle(origin int) Particle {
	eKinetic := 4. + rand.Float64()
	mu := rand.Float64()
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
	p.updateDynamicTimeMotionConstants(m)
	return p
}

func (p *Particle) getAxialEnergy() float64 {
	return p.eKinetic - p.trajectory.radialEnergy
}

func (p *Particle) getPotentialEnergy() float64 {
	return p.trajectory.totEnergy - p.eKinetic
}

func (trajectory *TrajectoryConstants) getStopPotential() float64 {
	return trajectory.radialEnergy - trajectory.totEnergy
}

func (trajectory *TrajectoryConstants) getTurnaroundX(m *Model) float64 {
	return m.LfromV(trajectory.getStopPotential())
}

type TrajectoryConstants struct {
	c1, c2, amplitude, phase, totEnergy, radialEnergy float64
}

func (p *Particle) updateDynamicTimeMotionConstants(m *Model) {
	var xSheathStart, eSheathStart float64
	if m.Parameters.CathodeFallLength < p.x {
		xSheathStart = m.Parameters.CathodeFallLength
		eSheathStart = p.getAxialEnergy() - m.VfromL(p.x) - m.VfromL(m.Parameters.CathodeFallLength)
	} else {
		xSheathStart = p.x
		eSheathStart = p.getAxialEnergy()
	}
	p.trajectory.c1 = xSheathStart - m.timeMotionP
	p.trajectory.c2 = math.Copysign(m.Parameters.CathodeFallLength*math.Sqrt(eSheathStart/m.timeMotionK), p.mu)
	p.trajectory.amplitude = math.Sqrt(p.trajectory.c1*p.trajectory.c1 + p.trajectory.c2*p.trajectory.c2)
	p.trajectory.phase = math.Atan(p.trajectory.c1 / p.trajectory.c2)
}

func (trajectory *TrajectoryConstants) SheathPositionAtTime(t float64, m *Model) (x float64) {
	return trajectory.amplitude*math.Sin(m.timeMotionW*t+trajectory.phase) + m.timeMotionP
}

func (trajectory *TrajectoryConstants) SheathTimeBetweenPositionsNoReversal(from, to float64, m *Model) (t float64) {
	fromTrig := (from - m.timeMotionP) / trajectory.amplitude
	toTrig := (to - m.timeMotionP) / trajectory.amplitude
	if math.Abs(fromTrig) > 1 {
		if math.Abs(fromTrig)-1 > 1e-5 {
			println("strange from")
		}
		fromTrig = math.Copysign(1., fromTrig)
	}
	if math.Abs(toTrig) > 1 {
		if math.Abs(toTrig)-1 > 1e-5 {
			println("strange to")
		}
		toTrig = math.Copysign(1., toTrig)
	}
	if to-from > 0 {
		t = (math.Asin(toTrig) - math.Asin(fromTrig)) / m.timeMotionW
	} else {
		t = -(math.Asin(toTrig) - math.Asin(fromTrig)) / m.timeMotionW
	}

	if t < 0 {
		println("t < 0")
	}
	return t
}

func (trajectory *TrajectoryConstants) ConstFieldTimeBetweenPositionsNoReversal(from, to float64, m *Model) (t float64) {
	vTo, vFrom := trajectory.axialVelocityAt(to, m), trajectory.axialVelocityAt(from, m)
	if to-from > 0 {
		t = (vTo - vFrom) / m.constFieldAcceleration
	} else {
		t = -(vTo - vFrom) / m.constFieldAcceleration
	}

	if t < 0 {
		println("t < 0")
	}
	return
}

func (trajectory *TrajectoryConstants) TimeBetweenPositionsNoReversal(from, to float64, m *Model) (t float64) {
	if from > m.Parameters.CathodeFallLength {
		if to > m.Parameters.CathodeFallLength {
			t = trajectory.ConstFieldTimeBetweenPositionsNoReversal(from, to, m)

			if math.IsNaN(t) {
				println("t is nan")
			}

			if t < 0 {
				println("t < 0")
			}
		} else {
			t = trajectory.ConstFieldTimeBetweenPositionsNoReversal(from, m.Parameters.CathodeFallLength, m) +
				trajectory.SheathTimeBetweenPositionsNoReversal(m.Parameters.CathodeFallLength, to, m)

			if math.IsNaN(t) {
				println("t is nan")
			}

			if t < 0 {
				println("t < 0")
			}
		}
	} else {
		if to > m.Parameters.CathodeFallLength {
			t = trajectory.SheathTimeBetweenPositionsNoReversal(from, m.Parameters.CathodeFallLength, m) +
				trajectory.ConstFieldTimeBetweenPositionsNoReversal(m.Parameters.CathodeFallLength, to, m)

			if math.IsNaN(t) {
				println("t is nan")
			}

			if t < 0 {
				println("t < 0")
			}
		} else {
			t = trajectory.SheathTimeBetweenPositionsNoReversal(from, to, m)

			if math.IsNaN(t) {
				println("t is nan")
			}

			if t < 0 {
				println("t < 0")
			}
		}
	}
	return t
}

func (trajectory *TrajectoryConstants) TFromX(p *Particle, x float64, future bool, m *Model) (t float64) {
	if future {
		if p.mu > 0 { //no reversal
			if p.x > x { // ---------x----------p.x-->>-----
				return math.NaN() //impossible
			} else { // ---------p.x-->>------x---------
				return trajectory.TimeBetweenPositionsNoReversal(p.x, x, m)
			}
		} else {
			turnaroundPoint := trajectory.getTurnaroundX(m)
			if p.x > x { // ---------x----------<<--p.x-----
				if turnaroundPoint < x { // ----*----x----------<<--p.x-----
					return trajectory.TimeBetweenPositionsNoReversal(p.x, x, m)
				} else { // ---------x----*-----<<--p.x-----
					return math.NaN()
				}
			} else { // ----*----<<--p.x---------x-----
				return trajectory.TimeBetweenPositionsNoReversal(p.x, turnaroundPoint, m) + trajectory.TimeBetweenPositionsNoReversal(turnaroundPoint, x, m)
			}
		}
	} else {
		if p.mu > 0 { // no reversal
			turnaroundPoint := trajectory.getTurnaroundX(m)
			if p.x > x { // ---------x---->>--p.x-----
				if turnaroundPoint < x { // ----*----x---->>--p.x-----
					return trajectory.TimeBetweenPositionsNoReversal(x, p.x, m)
				} else { // ----x---*---->>--p.x-----
					return math.NaN()
				}
			} else { // --*-->>--p.x--------x-----
				return trajectory.TimeBetweenPositionsNoReversal(x, turnaroundPoint, m) + trajectory.TimeBetweenPositionsNoReversal(turnaroundPoint, p.x, m)
			}
		} else {
			if p.x > x { // ----x------p.x--<<-----
				return math.NaN()
			} else { // ----p.x--<<----x---------
				return trajectory.TimeBetweenPositionsNoReversal(x, p.x, m)
			}
		}
	}
}

func (trajectory *TrajectoryConstants) SheathEnergyAtTime(t float64, m *Model) (e float64) {
	sp := trajectory.amplitude * math.Cos(m.timeMotionW*t+trajectory.phase) //C2*math.Cos(wt) - C1*math.Sin(wt)
	energyEV := m.timeMotionK * sp * sp * m.inverseCathodeFallLength * m.inverseCathodeFallLength
	return energyEV
}

func (trajectory *TrajectoryConstants) SheathTimeFromEnergy(e float64, m *Model) (t float64) {
	return (math.Acos(math.Sqrt(e*m.Parameters.CathodeFallLength*m.Parameters.CathodeFallLength/m.timeMotionK)/trajectory.amplitude) - trajectory.phase) / m.timeMotionW
}

func (trajectory *TrajectoryConstants) axialVelocityAt(x float64, m *Model) (v float64) {
	v = utils.EV2electronVelocity(trajectory.totEnergy + m.VfromL(x) - trajectory.radialEnergy)
	if math.IsNaN(v) {
		println("v is nan")
	}
	return
}

func (p *Particle) axialVelocityNow() (v float64) {
	return utils.EV2electronVelocity(p.getAxialEnergy())
}

func (p *Particle) setEnergy(eKinetic float64, s *Model, zeroChangeAcceptable bool, setX bool) {
	// p.MoveRadial(eKinetic, s, zeroChangeAcceptable)
	// r2 := p.y*p.y + p.z*p.z
	if eKinetic < p.trajectory.radialEnergy {
		fmt.Printf("stopPoint = %f\np: %v; %p\n", p.trajectory.getTurnaroundX(s), p, p)
		panic("eKinetic < p.trajectory.eStar")
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
	p.updateDynamicTimeMotionConstants(s)
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
	if p.trajectory.totEnergy > s.Va+s.Vc+5 {
		println("total energy greater than might ever be")
	}

	if s.Parameters.Volumetric {
		p.prevAxialEnergy = p.getAxialEnergy()
		p.prevX = p.x
		p.prevMuSign = p.mu
	}
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
	p.updateDynamicTimeMotionConstants(m)
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
