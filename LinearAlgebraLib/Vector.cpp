#include <iostream>
#include <cmath>
#include "Vector.h"
#include "constants.h"

using namespace std;

//DEFAULT CONSTRUCTOR
Vector::Vector()
{
	size = MAX;
	for (int i = 0; i < MAX; i++)
	{
		e[i] = 0.0;
		flag[i] = i;
	}
}

//OVERLOAD CONSTRUCTOR
Vector::Vector(int nSize)
{
	size = nSize;
	for (int i = 0; i < MAX; i++)
	{
		e[i] = 0.0;
		flag[i] = i;
	}
}


//DESTRUCTOR
Vector::~Vector()
{

}

//------------------------------------------------------------------------------
//PRINT
//------------------------------------------------------------------------------
void Vector::print() const
{
	for (int i = 0; i < size; i++)
		cout << "\t[" << i + 1 << "] = " << e[i] << "\n" << endl;
	cout << endl;
}

//------------------------------------------------------------------------------
//VECTOR ADDITION
//------------------------------------------------------------------------------
Vector Vector::Add(const Vector& other) const
{
	Vector R(size);
	for (int i = 0; i < size; i++)
		R.e[i] = e[i] + other.e[i];
	return R;
}

//------------------------------------------------------------------------------
//VECTOR SUBTRACTION
//------------------------------------------------------------------------------
Vector Vector::Sub(const Vector& other) const
{
	Vector R(size);
	for (int i = 0; i < size; i++)
		R.e[i] = e[i] - other.e[i];
	return R;
}

//------------------------------------------------------------------------------
//MULTIPLYING A VECTOR BY A SCALAR
//------------------------------------------------------------------------------
Vector Vector::scalarMulti(float factor) const
{
	Vector R(size);
	for (int i = 0; i < size; i++)
		R.e[i] = e[i] * factor;
	return R;
}

//------------------------------------------------------------------------------
//DIVIDE VECTOR BY SCALAR
//------------------------------------------------------------------------------
void Vector::divScalar(float den)
{
	for (int i = 0; i < size; i++)
		e[i] = e[i] / den;
}

//------------------------------------------------------------------------------
//DOT PRODUCT
//------------------------------------------------------------------------------
float Vector::Dot(const Vector& other) const
{
	float soma = 0.0;
	for (int i = 0; i < size; i++)
		soma += e[i] * other.e[i];
	return soma;
}

//------------------------------------------------------------------------------
//EUCLIDEAN NORM
//------------------------------------------------------------------------------
float Vector::Norm() const
{
	float soma = 0.0;
	for (int i = 0; i < size; i++)
		soma += pow(e[i], 2.0);
	return sqrt(soma);
}
