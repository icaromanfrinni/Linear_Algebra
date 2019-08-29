/*
 * BIBLIOTECA DE VETORES COM NÚMEROS COMPLEXOS
 */

#ifndef CVECTOR_H
#define CVECTOR_H

#include <ccomplex>
#include "constants.h"
#include "Vector.h"

using namespace std;

class cVector
{
public:
	//DEFAULT CONSTRUCTOR
	cVector();
	//OVERLOAD CONSTRUCTOR ('nSize' empty complex vector)
	cVector(int nSize);
	//OVERLOAD CONSTRUCTOR (from another complex vector)
	cVector(const cVector& other);
	//OVERLOAD CONSTRUCTOR (from a real vector)
	cVector(const Vector& realVector);
	//DESTRUCTOR
	~cVector();

	//--------------------------------
	//PRINT
	void print() const;
	//--------------------------------
	//MULTIPLYING A VECTOR BY A SCALAR
	cVector scalarMulti(complex<double> factor) const;
	//--------------------------------
	//DIVIDE A VECTOR BY A SCALAR
	void divScalar(complex<double> den);
	//--------------------------------
	//EUCLIDEAN NORM
	complex<double> Norm() const;
	//--------------------------------
	//DOT PRODUCT
	complex<double> Dot(const cVector& other) const;

	//--------------------------------
	//VARIABLES
	int size; //the number of rows
	complex<double> e[MAX];
	int flag[MAX];
};

#endif // CVECTOR_H
