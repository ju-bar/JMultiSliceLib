// FILE "V3D.h"
// HEADERFILE for "V3D.cpp"
// DECLARATION of class CV3D and related types and parameters
// AUTHOR Juri Barthel
// DATE 2014/12/31

#pragma once

#define V3D_CHOP_THR		1.0E-5

class CV3D
{
public:
	// construcors
	CV3D(void);
	CV3D(float x0, float y0, float z0);
	CV3D(const CV3D &v);
	// destructor
	~CV3D(void);

public:
	// member variables
	float x;
	float y;
	float z;

public:
	// operators
	//void operator=(CV3D v); // copies data
	void operator=(const CV3D &v); // copies data
	CV3D operator+(CV3D v); // adding
	CV3D operator+(const CV3D &v) const; // adding
	CV3D operator-(CV3D v); // substracting
	CV3D operator-(const CV3D &v) const; // substracting
	void operator+=(CV3D v); // add to
	void operator-=(CV3D v); // substract from
	CV3D operator*(float s); // scalar product
	CV3D operator/(float s); // scalar division
	const bool operator==(CV3D v) const; // compares data
	const bool operator!=(CV3D v) const; // compares data

	// functions
	float Item(int idx); // returns the vector item depending on index 0->x, 1->y, 2->z
	float Length(void);
	CV3D Normalized(void);
	float DotProduct(CV3D v); // dot product of "this" and "v"
	CV3D CrossProduct(CV3D v); // cross product of "this" and "v"
	bool Colinear(CV3D v); // checks colinearity with v
	CV3D MinCoords(CV3D v); // returns a vector with the minimum of the coordinates of this and v
	CV3D MaxCoords(CV3D v); // returns a vector with the maximum of the coordinates of this and v
	CV3D Scale(CV3D scale); // returns a scaled version of this
	CV3D Floor(void); // returns a vector with all coordinates set to the floor of the coordinates of this
	CV3D Ceil(void); // returns a vector with all coordinates set to the ceiling of the coordinates of this
	CV3D Round(void); // returns a vector with all coordinates set to the next integer of the coordinates of this
	void Chop(float fRelChopThr=V3D_CHOP_THR); // Chop small numbers when their contribution in length is small
	CV3D Modulo(CV3D vModulo); // returns a vector with coordinates modulo another vector
};

