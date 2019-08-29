#include <iostream>
#include <cmath>
#include "cVector.h"

using namespace std;

//DEFAULT CONSTRUCTOR
cVector::cVector()
{
	size = MAX;
	for (int i = 0; i < MAX; i++)
	{
		e[i] = 0.0;
		flag[i] = i;
	}
}

//OVERLOAD CONSTRUCTOR ('nSize' empty complex vector)
cVector::cVector(int nSize)
{
	size = nSize;
	for (int i = 0; i < MAX; i++)
	{
		e[i] = 0.0;
		flag[i] = i;
	}
}

//OVERLOAD CONSTRUCTOR (from another complex vector)
cVector::cVector(const cVector& other)
{
	size = other.size;
	for (int i = 0; i < MAX; i++)
	{
		e[i] = other.e[i];
		flag[i] = other.flag[i];
	}
}

//OVERLOAD CONSTRUCTOR (from a real vector)
cVector::cVector(const Vector& realVector)
{
	size = realVector.size;
	for (int i = 0; i < MAX; i++)
	{
		e[i] = realVector.e[i];
		flag[i] = realVector.flag[i];
	}
}

//DESTRUCTOR
cVector::~cVector()
{
}

//------------------------------------------------------------------------------
//PRINT
//------------------------------------------------------------------------------
void cVector::print() const
{
	for (int i = 0; i < size; i++)
		cout << "\t[" << i + 1 << "] = " << e[i].real() << " " << showpos << e[i].imag() << "i" << noshowpos << "\n";
	cout << endl;
}

//------------------------------------------------------------------------------
//MULTIPLYING A VECTOR BY A SCALAR
//------------------------------------------------------------------------------
cVector cVector::scalarMulti(complex<double> factor) const
{
	cVector R(size);
	for (int i = 0; i < size; i++)
		R.e[i] = e[i] * factor;
	return R;
}

//------------------------------------------------------------------------------
//DIVIDE VECTOR BY SCALAR
//------------------------------------------------------------------------------
void cVector::divScalar(complex<double> den)
{
	for (int i = 0; i < size; i++)
		e[i] = e[i] / den;
}

//------------------------------------------------------------------------------
//DOT PRODUCT
//------------------------------------------------------------------------------
complex<double> cVector::Dot(const cVector& other) const
{
	complex<double> soma = 0.0;
	for (int i = 0; i < size; i++)
		soma += e[i] * conj(other.e[i]);
	return soma;
}

//------------------------------------------------------------------------------
//EUCLIDEAN NORM
//------------------------------------------------------------------------------
complex<double> cVector::Norm() const
{
	complex<double> soma = 0.0;
	for (int i = 0; i < size; i++)
		soma += norm(e[i]);
	return sqrt(soma);
}
