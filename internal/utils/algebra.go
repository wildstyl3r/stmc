package utils

import "fmt"

type Vec2D struct {
	X, Y float64
}

func (v Vec2D) String() string {
	return fmt.Sprintf("[%f,%f]", v.X, v.Y)
}

func (l Vec2D) Add(r Vec2D) Vec2D {
	return Vec2D{l.X + r.X, l.Y + r.Y}
}
func (l Vec2D) Sub(r Vec2D) Vec2D {
	return Vec2D{l.X - r.X, l.Y - r.Y}
}
func (l Vec2D) Mul(r float64) Vec2D {
	return Vec2D{l.X * r, l.Y * r}
}
func (l Vec2D) Div(r float64) Vec2D {
	return Vec2D{l.X / r, l.Y / r}
}

func (l Vec2D) PointMul(r Vec2D) Vec2D {
	return Vec2D{l.X * r.X, l.Y * r.Y}
}

func (v Vec2D) InsideOfRectangle(lowerLeft, upperRight Vec2D) bool {
	return lowerLeft.X < v.X && v.X < upperRight.X && lowerLeft.Y < v.Y && v.Y < upperRight.Y
}

func Vecs2Dto2slices(vs []Vec2D) (sx, sy []float64) {
	sx, sy = make([]float64, len(vs)), make([]float64, len(vs))
	for i := range vs {
		sx[i] = vs[i].X
		sy[i] = vs[i].Y
	}
	return
}

func SlicesToVec2D(xs, ys []float64) (vs []Vec2D) {
	vs = make([]Vec2D, min(len(xs), len(ys)))
	for i := range min(len(xs), len(ys)) {
		vs[i].X = xs[i]
		vs[i].Y = ys[i]
	}
	return
}

type Mat2D struct {
	a, b float64
	c, d float64
}

func (m *Mat2D) Det() float64 {
	return m.a*m.d - m.b*m.c
}

func (m *Mat2D) Inverse() Mat2D {
	invDet := 1 / m.Det()
	return Mat2D{
		m.d * invDet, -m.b * invDet,
		-m.c * invDet, m.a * invDet}
}

func (l Mat2D) MulVec(r Vec2D) Vec2D {
	return Vec2D{
		l.a*r.X + l.b*r.Y,
		l.c*r.X + l.d*r.Y}
}

// type Mat3D struct {
// 	a, b, c float64
// 	d, e, f float64
// 	g, h, i float64
// }

// func (m *Mat3D) Det() float64 {
// 	return m.a*m.e*m.i + m.b*m.f*m.g + m.d*m.h*m.c - m.c*m.e*m.g - m.b*m.d*m.i - m.f*m.h*m.a
// }

// func (m *Mat3D) Inverse() Mat3D {
// 	invDet := 1 / m.Det()
// 	return Mat3D{
// 		m.d * invDet, -m.b * invDet,
// 		-m.c * invDet, m.a * invDet}
// }

// func (l Mat3D) MulVec(r Vec3D) Vec3D {
// 	return Vec3D{
// 		l.a*r.X + l.b*r.Y,
// 		l.c*r.X + l.d*r.Y}
// }
