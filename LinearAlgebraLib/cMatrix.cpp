#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <cmath>
#include "cMatrix.h"

using namespace std;

//DEFAULT CONSTRUCTOR
cMatrix::cMatrix()
{
	sizeR = 0;
	sizeC = 0;
	for (int i = 0; i < MAX; i++)
		for (int j = 0; j < MAX; j++)
			e[i][j] = 0.0;
}

//OVERLOAD CONSTRUCTOR
cMatrix::cMatrix(int nR, int nC)
{
	sizeR = nR;
	sizeC = nC;
	for (int i = 0; i < MAX; i++)
		for (int j = 0; j < MAX; j++)
			e[i][j] = 0.0;
}


//DESTRUCTOR
cMatrix::~cMatrix()
{

}


//------------------------------------------------------------------------------
//PRINT
//------------------------------------------------------------------------------
void cMatrix::print() const
{
	for (int i = 0; i < sizeR; i++)
	{
		for (int j = 0; j < sizeC; j++)
			cout << "\t" << e[i][j].real() << " " << showpos << e[i][j].imag() << "i" << noshowpos;
		cout << "\n" << endl;
	}
}

//------------------------------------------------------------------------------
//SCALAR MULTIPLICATION
//------------------------------------------------------------------------------
cMatrix cMatrix::scalarMulti(complex<double> factor) const
{
	cMatrix R(sizeR, sizeC);
	for (int i = 0; i < sizeR; i++)
		for (int j = 0; j < sizeC; j++)
			R.e[i][j] = e[i][j] * factor;
	return R;
}

//------------------------------------------------------------------------------
//MULTIPLYING A MATRIX BY A VECTOR
//------------------------------------------------------------------------------
cVector cMatrix::vectorMulti(const cVector &v) const
{
	if (sizeC != v.size)
	{
		cerr << "\n!!! INNER MATRIX DIMENSIONS MUST AGREE !!!\n" << endl;
		exit(EXIT_FAILURE);
	}
	cVector R(sizeR);
	for (int j = 0; j < sizeC; j++)
		for (int i = 0; i < sizeR; i++)
			R.e[i] += e[i][j] * v.e[j];
	return R;
}

//------------------------------------------------------------------------------
//MULTIPLYING A MATRIX BY A MATRIX
//------------------------------------------------------------------------------
cMatrix cMatrix::matrixMulti(const cMatrix &other) const
{
	if (sizeC != other.sizeR)
	{
		cerr << "\n!!! INNER MATRIX DIMENSIONS MUST AGREE !!!\n" << endl;
		exit(EXIT_FAILURE);
	}
	cMatrix R(sizeR, other.sizeC);
	for (int j = 0; j < other.sizeC; j++)
		for (int k = 0; k < sizeR; k++)
			for (int i = 0; i < sizeR; i++)
				R.e[i][j] += e[i][k] * other.e[k][j];
	return R;
}

//------------------------------------------------------------------------------
//AUGMENTED MATRIX (MATRIX WITH MATRIX)
//------------------------------------------------------------------------------
cMatrix cMatrix::augMatrix(const cMatrix& other) const
{
	cMatrix aM(sizeR, sizeC + other.sizeC);
	for (int i = 0; i < sizeR; i++)
	{
		for (int j = 0; j < sizeC; j++)
			aM.e[i][j] = e[i][j];
		for (int j = sizeC; j < aM.sizeC; j++)
			aM.e[i][j] = other.e[i][j - sizeC];
	}
	return aM;
}

//------------------------------------------------------------------------------
//IDENTITY MATRIX
//------------------------------------------------------------------------------
cMatrix IdentC(int size)
{
	cMatrix I(size, size);
	for (int i = 0; i < size; i++)
		I.e[i][i] = 1;
	return I;
}
