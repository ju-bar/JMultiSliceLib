// FILE "M3D.cpp"
// IMPLEMENTATION of class CM3D
// AUTHOR Juri Barthel
// DATE 2016/01/18

#include "M3D.h"
#include <math.h>

CM3D::CM3D(void)
{
	Identity();
}

CM3D::CM3D(CV3D x0, CV3D y0, CV3D z0)
{
	x = x0;
	y = y0;
	z = z0;
}

CM3D::CM3D(const CM3D &m)
{
	x = m.x;
	y = m.y;
	z = m.z;
}

CM3D::~CM3D(void)
{
}

void CM3D::operator=(const CM3D &m)
{
	x = m.x;
	y = m.y;
	z = m.z;
}

CM3D CM3D::operator*(float s)
{
	return CM3D(x*s,y*s,z*s);
}

CM3D CM3D::operator/(float s)
{
	if (s==0.0f)
		return *this;
	return CM3D(x/s,y/s,z/s);
}

const bool CM3D::operator==(CM3D m) const
{
	return ((x == m.x) && (y == m.y) && (z == m.z));
}

const bool CM3D::operator!=(CM3D m) const
{
	return ((x != m.x) || (y != m.y) || (y != m.y));
}

void CM3D::Identity(void)
{
	x = CV3D(1.0f,0.0f,0.0f);
	y = CV3D(0.0f,1.0f,0.0f);
	z = CV3D(0.0f,0.0f,1.0f);
}

float CM3D::Determinant(void)
{
	return (x.x*y.y-x.y*y.x)*z.z+(x.y*y.z-y.y*x.z)*z.x+(x.z*y.x-x.x*y.z)*z.y;
}

float CM3D::Trace(void)
{
	return x.x+y.y+z.z;
}

void CM3D::Normalize(void)
{
	float d = Determinant();
	if (d!=0.0f) {
		x = x/d;
		y = y/d;
		z = z/d;
	}
}

CM3D CM3D::Normalized(void)
{
	CM3D m(*this);
	m.Normalize();
	return m;
}

void CM3D::Transpose(void)
{
	float ftmp = 0.0f;
	ftmp = x.y; x.y = y.x; y.x = ftmp;
	ftmp = x.z; x.z = z.x; z.x = ftmp;
	ftmp = y.z; y.z = z.y; z.y = ftmp;
}

CM3D CM3D::Transposed(void)
{
	CM3D m(*this);
	m.Transpose();
	return m;
}

CV3D CM3D::Item(int idx)
{
	if (idx<=0) return x;
	if (idx==1) return y;
	if (idx>=2) return z;
	return CV3D();
}

float CM3D::Item(int idx, int jdx)
{
	return Item(idx).Item(jdx);
}

void CM3D::Invert(void)
{
	CM3D m;
	m = (*this).Inverse();
	*this = m;
}

CM3D CM3D::Inverse(void)
{
	CM3D m;
	float det = Determinant();
	if (det!=0.f) {
		float fidet = 1.f/det;
		m.x.x = (y.y*z.z - z.y*y.z)*fidet;
		m.x.y = (z.y*x.z - x.y*z.z)*fidet;
		m.x.z = (x.y*y.z - y.y*x.z)*fidet;
		m.y.x = (y.z*z.x - z.z*y.x)*fidet;
		m.y.y = (z.z*x.x - x.z*z.x)*fidet;
		m.y.z = (x.z*y.x - y.z*x.x)*fidet;
		m.z.x = (y.x*z.y - z.x*y.y)*fidet;
		m.z.y = (z.x*x.y - x.x*z.y)*fidet;
		m.z.z = (x.x*y.y - y.x*x.y)*fidet;
	}
	return m;
}

CV3D CM3D::VProduct(CV3D v)
{
	return CV3D( x.DotProduct(v), y.DotProduct(v), z.DotProduct(v) );
}

CM3D CM3D::MProduct(CM3D m)
{
	CM3D tm = m;
	tm.Transpose();
	CM3D mm;
	mm.x.x = x.DotProduct(tm.x);
	mm.x.y = x.DotProduct(tm.y);
	mm.x.z = x.DotProduct(tm.z);
	mm.y.x = y.DotProduct(tm.x);
	mm.y.y = y.DotProduct(tm.y);
	mm.y.z = y.DotProduct(tm.z);
	mm.z.x = z.DotProduct(tm.x);
	mm.z.y = z.DotProduct(tm.y);
	mm.z.z = z.DotProduct(tm.z);
	return mm;
}

void CM3D::Chop(float fRelChopThr)
{
	x.Chop(fRelChopThr);
	y.Chop(fRelChopThr);
	z.Chop(fRelChopThr);
}