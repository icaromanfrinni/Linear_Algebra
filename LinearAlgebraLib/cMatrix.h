/*
 * BIBLIOTECA DE OPERAÇÕES MATRICIAIS COM NÚMEROS COMPLEXOS
 */

#ifndef CMATRIX_H
#define CMATRIX_H

#include "constants.h"
#include "cVector.h"

class cMatrix
{
public:
	//DEFAULT CONSTRUCTOR
	cMatrix();
	//OVERLOAD CONSTRUCTOR
	cMatrix(int nR, int nC);
	//DESTRUCTOR
	~cMatrix();

	//--------------------------------
	//PRINT
	void print() const;
	//--------------------------------
	//SCALAR MULTIPLICATION
	cMatrix scalarMulti(complex<double> factor) const;
	//--------------------------------
	//MULTIPLYING A MATRIX BY A VECTOR
	cVector vectorMulti(const cVector& v) const;
	//--------------------------------
	//MULTIPLYING A MATRIX BY A MATRIX
	cMatrix matrixMulti(const cMatrix& other) const;
	//--------------------------------
	//AUGMENTED MATRIX (MATRIX WITH MATRIX)
	cMatrix augMatrix(const cMatrix& other) const;

	//--------------------------------
	//VARIABLES
	int sizeR; //the number of rows
	int sizeC; //the number of columns
	complex<double> e[MAX][MAX];
};

//====================================
//IDENTITY MATRIX
cMatrix IdentC(int size);

#endif // CMATRIX_H
