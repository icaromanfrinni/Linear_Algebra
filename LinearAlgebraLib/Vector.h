/*
 * BIBLIOTECA DE VETORES
 */

#ifndef VECTOR_H
#define VECTOR_H

#include "constants.h"

class Vector
{
public:
	//DEFAULT CONSTRUCTOR
	Vector();
	//OVERLOAD CONSTRUCTOR
	Vector(int nSize);
	//DESTRUCTOR
	~Vector();

	//--------------------------------
	//PRINT
	void print() const;
	//--------------------------------
	//VECTOR ADDITION
	Vector Add(const Vector& other) const;
	//--------------------------------
	//VECTOR SUBTRACTION
	Vector Sub(const Vector& other) const;
	//--------------------------------
	//MULTIPLYING A VECTOR BY A SCALAR
	Vector scalarMulti(float factor) const;
	/*
	//--------------------------------
	//MULTIPLYING A VECTOR BY A VECTOR
	Vector vectorMulti(const Vector& other) const;
	*/
	//--------------------------------
	//DIVIDE A VECTOR BY A SCALAR
	void divScalar(float den);
	//--------------------------------
	//EUCLIDEAN NORM
	float Norm() const;
	//--------------------------------
	//DOT PRODUCT
	float Dot(const Vector& other) const;

	//--------------------------------
	//VARIABLES
	int size; //the number of rows
	double e[MAX];
	int flag[MAX];
};

#endif // VECTOR_H
