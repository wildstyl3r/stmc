package utils

import (
	"fmt"
	"math"
)

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

func GetLinearCoefficients(x1, y1, x2, y2 float64) (a, b float64) {
	a = (y2 - y1) / (x2 - x1)
	b = y1 - x1*a
	return a, b
}

type Vec3D struct {
	X, Y, Z float64
}

func (a *Vec3D) Add(b *Vec3D) Vec3D {
	return Vec3D{
		a.X + b.X,
		a.Y + b.Y,
		a.Z + b.Z,
	}
}

func (a *Vec3D) Sub(b *Vec3D) Vec3D {
	return Vec3D{
		a.X - b.X,
		a.Y - b.Y,
		a.Z - b.Z,
	}
}

func (a *Vec3D) Neg() Vec3D {
	return Vec3D{
		-a.X,
		-a.Y,
		-a.Z,
	}
}

func (a *Vec3D) Dot(b *Vec3D) float64 {
	return a.X*b.X + a.Y*b.Y + a.Z*b.Z
}

func (a *Vec3D) Norm() float64 {
	return math.Sqrt(a.Dot(a))
}

func (a *Vec3D) Norm2() float64 {
	return a.Dot(a)
}

func (a *Vec3D) Scale(scale float64) Vec3D {
	return Vec3D{
		a.X * scale,
		a.Y * scale,
		a.Z * scale,
	}
}

func (a *Vec3D) Unit() Vec3D {
	n := a.Norm()
	return Vec3D{
		a.X / n,
		a.Y / n,
		a.Z / n,
	}
}

func (a *Vec3D) Cross(b *Vec3D) Vec3D {
	return Vec3D{
		a.Y*b.Z - a.Z*b.Y,
		a.Z*b.X - a.X*b.Z,
		a.X*b.Y - a.Y*b.X,
	}
}

// http://jcgt.org/published/0006/01/01/
func (v *Vec3D) ModifiedFrisvadBasis() (b1 Vec3D, b2 Vec3D) {
	s := math.Copysign(1, v.Z)
	a := -1. / (s + v.Z)
	b := v.X * v.Y * a
	b1 = Vec3D{
		1 + s*v.X*v.X*a,
		s * b,
		-s * v.X,
	}
	b2 = Vec3D{
		b,
		s + v.Y*v.Y*a,
		-v.Y,
	}
	return
}

func RodriguesRotation(v, axis *Vec3D, angle float64) (result Vec3D) {
	k := axis.Unit()
	sin, cos := math.Sincos(angle)

	vrej := k.Cross(v)
	vrej = vrej.Scale(sin)

	vproj := k.Scale(k.Dot(v) * (1 - cos))

	vscaled := v.Scale(cos)
	result = vscaled.Add(&vrej)
	result = result.Add(&vproj)
	return result
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
