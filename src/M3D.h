// FILE "M3D.h"
// HEADERFILE for "M3D.cpp"
// DECLARATION of class CM3D and related types and parameters
// AUTHOR Juri Barthel
// DATE 2016/01/19

#pragma once

#include "V3D.h"

class CM3D
{
public:
	// construcors
	CM3D(void);
	CM3D(CV3D x0, CV3D y0, CV3D z0);
	CM3D(const CM3D &m);
	// destructor
	~CM3D(void);

public:
	// member variables
	CV3D x;
	CV3D y;
	CV3D z;

public:
	// operators
	void operator=(const CM3D &m); // copies data
	CM3D operator*(float s); // scalar product
	CM3D operator/(float s); // scalar division
	const bool operator==(CM3D m) const; // compares data
	const bool operator!=(CM3D m) const; // compares data

	// functions
	CV3D Item(int idx); // returns row item idx
	float Item(int idx, int jdx); // returns column item jdx in row item idx
	void Identity(void); // Initializes "this" as identity matrix
	float Determinant(void); // Returns the determinant of "this"
	float Trace(void); // Returns the trace of "this"
	void Normalize(void); // Normalizes "this" (determinant = 1)
	CM3D Normalized(void); // Returns a normalized version of "this" (determinant = 1)
	void Transpose(void); // Transposes "this"
	CM3D Transposed(void); // Returns a transposed version of "this"
	void Invert(void); // Inverts "this"
	CM3D Inverse(void); // Returns the inverse of "this"
	CV3D VProduct(CV3D v); // Product of "this" and vector "v"
	CM3D MProduct(CM3D m); // product of "this" and matrix "m"
	void Chop(float fRelChopThr=V3D_CHOP_THR); // Chop small numbers when their contribution in length is small
};

