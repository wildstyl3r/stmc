package model

import (
	"math"

	"github.com/wildstyl3r/stmc/internal/utils"
)

type TrajectoryConstants struct {
	radialEnergy float64 //[eV]
	totEnergy    float64 //aka Vcap
}

func (trj *TrajectoryConstants) TimeFromReversal(x float64, m *Model) (t float64) {
	cos := (x - m.motionK) / (trj.getTurnaroundX(m) - m.motionK)
	if math.Abs(cos) > 1 {
		cos = math.Copysign(1., cos)
	}
	return math.Acos(cos) / m.motionW
}

func (trj *TrajectoryConstants) PositionAtTimeFromReversal(t float64, m *Model) (x float64) {
	return (trj.getTurnaroundX(m)-m.motionK)*math.Cos(m.motionW*t) + m.motionK
}

func (trj *TrajectoryConstants) getStopPotential() float64 {
	return trj.radialEnergy - trj.totEnergy
}

func (trj *TrajectoryConstants) getTurnaroundX(m *Model) float64 {
	return m.LfromV(trj.getStopPotential())
}
func (trj *TrajectoryConstants) getAxialEnergyAt(position float64, m *Model) float64 {
	return trj.totEnergy + m.VfromL(position) - trj.radialEnergy
}

func (trj *TrajectoryConstants) SheathPositionAfterTime(startX, t float64, m *Model, prograde bool) (x float64, endPrograde bool) {
	tAsRev := trj.TimeFromReversal(startX, m)
	if prograde {
		return trj.PositionAtTimeFromReversal(t+tAsRev, m), true
	} else {
		if t > tAsRev {
			return trj.PositionAtTimeFromReversal(t-tAsRev, m), true
		} else {
			return trj.PositionAtTimeFromReversal(tAsRev-t, m), false
		}
	}
}

func (trj *TrajectoryConstants) ConstFieldPositionAfterTime(startX, t float64, m *Model, prograde bool) (x float64, endPrograde bool) {
	if prograde {
		return startX + trj.axialVelocityAt(startX, m)*t + 0.5*m.constFieldAcceleration*t*t, true
	} else {
		return startX - trj.axialVelocityAt(startX, m)*t + 0.5*m.constFieldAcceleration*t*t, -trj.axialVelocityAt(startX, m)+m.constFieldAcceleration*t > 0
	}
}

func (trj *TrajectoryConstants) SheathTimeToPosition(startX, endX float64, m *Model, prograde bool) (t float64) {
	if prograde {
		return trj.TimeFromReversal(endX, m) - trj.TimeFromReversal(startX, m)
	} else {
		return trj.TimeFromReversal(startX, m) - trj.TimeFromReversal(endX, m)
	}
}

func (trj *TrajectoryConstants) TurnaroundTime(startX float64, m *Model) (t float64) {
	if startX > m.Parameters.CathodeFallLength {
		return trj.TimeFromReversal(m.Parameters.CathodeFallLength, m) + trj.ConstFieldTimeToPosition(startX, m.Parameters.CathodeFallLength, m, false)
	} else {
		return trj.TimeFromReversal(startX, m)
	}
}

func (trj *TrajectoryConstants) ConstFieldTimeToPosition(startX, endX float64, m *Model, prograde bool) (t float64) {
	vTo, vFrom := trj.axialVelocityAt(endX, m), trj.axialVelocityAt(startX, m)
	if prograde {
		t = (vTo - vFrom) / m.constFieldAcceleration
	} else {
		t = (vFrom - vTo) / m.constFieldAcceleration
	}

	if t < 0 {
		println("t < 0")
	}
	return
}

func (trj *TrajectoryConstants) PositionAfterTime(startX, t float64, m *Model, prograde bool) (x float64, endPrograde bool) {
	if !prograde {
		turnaroundX := trj.getTurnaroundX(m)
		var turnaroundTime float64
		if turnaroundX < m.Parameters.CathodeFallLength {
			turnaroundTime = trj.TurnaroundTime(startX, m)
			if turnaroundTime > t {
				if startX < m.Parameters.CathodeFallLength {
					x, endPrograde = trj.SheathPositionAfterTime(startX, t, m, false)
					if math.IsNaN(x) {
						println("x is nan")
					}
					if x < turnaroundX {
						x = turnaroundX
					}
					return
				} else {
					tToSheathBoundary := trj.TimeBetweenPositionsNoReversal(startX, m.Parameters.CathodeFallLength, m, false)
					if tToSheathBoundary > t {
						x, endPrograde = trj.ConstFieldPositionAfterTime(startX, t, m, false)
						if math.IsNaN(x) {
							println("x is nan")
						}
						if x < turnaroundX {
							x = turnaroundX
						}
						return
					} else {
						x, endPrograde = trj.SheathPositionAfterTime(m.Parameters.CathodeFallLength, t-tToSheathBoundary, m, false)
						if math.IsNaN(x) {
							println("x is nan")
						}
						if x < turnaroundX {
							x = turnaroundX
						}
						return
					}
				}
			} else {
				t -= turnaroundTime
				prograde = true
				startX = turnaroundX
			}
		} else {
			x, endPrograde = trj.ConstFieldPositionAfterTime(startX, t, m, false)
			if math.IsNaN(x) {
				println("x is nan")
			}
			if x < turnaroundX {
				x = turnaroundX
			}
			return
		}
	}
	if startX > m.Parameters.CathodeFallLength {
		x, endPrograde = trj.ConstFieldPositionAfterTime(startX, t, m, true)
		if math.IsNaN(x) {
			println("x is nan")
		}
		if x < trj.getTurnaroundX(m) {
			x = trj.getTurnaroundX(m)
		}
		return
	} else {
		tToD := trj.TimeFromReversal(m.Parameters.CathodeFallLength, m)
		timeDelta := t - tToD
		if timeDelta > 0 {
			x, endPrograde = trj.ConstFieldPositionAfterTime(m.Parameters.CathodeFallLength, timeDelta, m, true)
			if math.IsNaN(x) {
				println("x is nan")
			}
			if x < trj.getTurnaroundX(m) {
				x = trj.getTurnaroundX(m)
			}
			return
		} else {
			x, endPrograde = trj.SheathPositionAfterTime(startX, t, m, true)
			if math.IsNaN(x) {
				println("x is nan")
			}
			if x < trj.getTurnaroundX(m) {
				x = trj.getTurnaroundX(m)
			}
			return
		}
	}
}

func (trj *TrajectoryConstants) TimeBetweenPositionsNoReversal(startX, endX float64, m *Model, prograde bool) (t float64) {
	if startX < m.Parameters.CathodeFallLength {
		if endX < m.Parameters.CathodeFallLength {
			t = trj.SheathTimeToPosition(startX, endX, m, prograde)
		} else {
			t = trj.SheathTimeToPosition(startX, m.Parameters.CathodeFallLength, m, prograde) + trj.ConstFieldTimeToPosition(m.Parameters.CathodeFallLength, endX, m, prograde)
		}
	} else {
		if endX < m.Parameters.CathodeFallLength {
			t = trj.ConstFieldTimeToPosition(startX, m.Parameters.CathodeFallLength, m, prograde) + trj.SheathTimeToPosition(m.Parameters.CathodeFallLength, endX, m, prograde)
		} else {
			trj.ConstFieldTimeToPosition(startX, endX, m, prograde)
		}
	}
	return t
}

func (trj *TrajectoryConstants) axialVelocityAt(x float64, m *Model) (v float64) {
	v = utils.EV2electronVelocity(max(0, trj.totEnergy+m.VfromL(x)-trj.radialEnergy))
	if math.IsNaN(v) {
		println("v is nan")
	}
	return
}
