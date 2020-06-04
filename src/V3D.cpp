// FILE "V3D.cpp"
// IMPLEMENTATION of class CV3D
// AUTHOR Juri Barthel
// DATE 2014/12/31

#include "V3D.h"
#include <math.h>
#include <algorithm>

CV3D::CV3D(void)
{
	x = 0.0f;
	y = 0.0f;
	z = 0.0f;
}

CV3D::CV3D(float x0, float y0, float z0)
{
	x = x0;
	y = y0;
	z = z0;
}

CV3D::CV3D(const CV3D &v)
{
	x = v.x;
	y = v.y;
	z = v.z;
}

CV3D::~CV3D(void)
{
}

void CV3D::operator=(const CV3D &v)
{
	x = v.x;
	y = v.y;
	z = v.z;
}

CV3D CV3D::operator+(CV3D v)
{
	return CV3D(x+v.x,y+v.y,z+v.z);
}

CV3D CV3D::operator+(const CV3D &v) const
{
	return CV3D(x+v.x,y+v.y,z+v.z);
}

CV3D CV3D::operator-(CV3D v)
{
	return CV3D(x-v.x,y-v.y,z-v.z);
}

CV3D CV3D::operator-(const CV3D &v) const
{
	return CV3D(x-v.x,y-v.y,z-v.z);
}


void CV3D::operator+=(CV3D v)
{
	x += v.x;
	y += v.y;
	z += v.z;
	return;
}

void CV3D::operator-=(CV3D v)
{
	x -= v.x;
	y -= v.y;
	z -= v.z;
	return;
}

CV3D CV3D::operator*(float s)
{
	return CV3D(x*s,y*s,z*s);
}

CV3D CV3D::operator/(float s)
{
	if (s==0.0f)
		return *this;
	return CV3D(x/s,y/s,z/s);
}

const bool CV3D::operator==(CV3D v) const
{
	return ((x == v.x) && (y == v.y) && (z == v.z));
}

const bool CV3D::operator!=(CV3D v) const
{
	return ((x != v.x) || (y != v.y) || (z != v.z));
}

float CV3D::Item(int idx)
{
	if (idx<=0) return x;
	if (idx==1) return y;
	if (idx>=2) return z;
	return 0.0f;
}

float CV3D::Length(void)
{
	return sqrtf(x*x+y*y+z*z);
}

CV3D CV3D::Normalized(void)
{
	float l = Length();
	if (l==0.0f) 
		return *this;
	return CV3D(*this)/l;
}

float CV3D::DotProduct(CV3D v)
{
	return x*v.x+y*v.y+z*v.z;
}

CV3D CV3D::CrossProduct(CV3D v)
{
	return CV3D(y*v.z-z*v.y,z*v.x-x*v.z,x*v.y-y*v.x);
}

bool CV3D::Colinear(CV3D v)
{
	CV3D v1 = this->Normalized();
	CV3D v2 = v.Normalized();
	return (1.0f==v1.DotProduct(v2));
}

CV3D CV3D::MinCoords(CV3D v)
{
	CV3D result;
	result.x = std::min(x,v.x);
	result.y = std::min(y,v.y);
	result.z = std::min(z,v.z);
	return result;
}

CV3D CV3D::MaxCoords(CV3D v)
{
	CV3D result;
	result.x = std::max(x,v.x);
	result.y = std::max(y,v.y);
	result.z = std::max(z,v.z);
	return result;
}

CV3D CV3D::Scale(CV3D scale)
{
	return CV3D(x*scale.x,y*scale.y,z*scale.z);
}

CV3D CV3D::Floor(void)
{
	return CV3D(floorf(x),floorf(y),floorf(z));
}

CV3D CV3D::Ceil(void)
{
	return CV3D(ceilf(x),ceilf(y),ceilf(z));
}

CV3D CV3D::Round(void)
{
	return CV3D(floorf(x+0.5f),floorf(y+0.5f),floorf(z+0.5f));
}

void CV3D::Chop(float fRelChopThr)
{
	float vlen = Length();
	if (vlen>0.0f) {
		if ( fabsf(x/vlen) < fRelChopThr ) x = 0.0f;
		if ( fabsf(y/vlen) < fRelChopThr ) y = 0.0f;
		if ( fabsf(z/vlen) < fRelChopThr ) z = 0.0f;
	}
}

CV3D CV3D::Modulo(CV3D vModulo)
{
	CV3D vResult = CV3D(0.f,0.f,0.f);
	float a, p;
	a = x; p = vModulo.x;
	if (p!=0.f) vResult.x = a - floorf( a / p) * p;
	a = y; p = vModulo.y;
	if (p!=0.f) vResult.y = a - floorf( a / p) * p;
	a = z; p = vModulo.z;
	if (p!=0.f) vResult.z = a - floorf( a / p) * p;
	return vResult;
}