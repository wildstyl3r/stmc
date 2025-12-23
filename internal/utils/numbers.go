package utils

import "math"

func IntAbs(a int) int {
	if a < 0 {
		return -a
	} else {
		return a
	}

}

func KahanSum(sum, piece, compensation *float64) {
	y := *piece - *compensation
	temp := *sum + y
	*compensation = (temp - *sum) - y
	*sum = temp
}

type KahanSummable struct {
	val, compensation float64
}

func (ks *KahanSummable) Add(piece float64) {
	// y := piece - ks.compensation
	// temp := ks.val + y
	// ks.compensation = (temp - ks.val) - y
	// ks.val = temp

	if math.IsNaN(piece) {
		println("piece is nan")
		return
	}

	t := ks.val + piece
	if math.Abs(ks.val) >= math.Abs(piece) {
		ks.compensation += (ks.val - t) + piece
	} else {
		ks.compensation += (piece - t) + ks.val
	}
	ks.val = t
}

func (ks *KahanSummable) Val() float64 {
	return ks.val + ks.compensation
}

func NewKahanSummable(from float64) KahanSummable {
	return KahanSummable{
		val: from,
	}
}
