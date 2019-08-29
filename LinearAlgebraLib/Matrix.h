/*
 * BIBLIOTECA DE OPERAÃ‡Ã•ES MATRICIAIS
 */

#ifndef MATRIX_H
#define MATRIX_H

#include "constants.h"
#include "Vector.h"
#include "cMatrix.h"

class Matrix
{
public:
	//DEFAULT CONSTRUCTOR
	Matrix();
	//OVERLOAD CONSTRUCTOR
	Matrix(int nR, int nC);
	//DESTRUCTOR
	~Matrix();

	//--------------------------------
	//PRINT
	void print() const;
	//--------------------------------
	//MATRIX ADDITION
	Matrix matrixAdd(const Matrix& other) const;
	//--------------------------------
	//Complex MATRIX ADDITION
	cMatrix matrixAdd(const cMatrix& other) const;
	//--------------------------------
	//SCALAR MULTIPLICATION
	Matrix scalarMulti(float factor) const;
	//--------------------------------
	//MULTIPLYING A MATRIX BY A VECTOR
	Vector vectorMulti(const Vector& v) const;
	//--------------------------------
	//MULTIPLYING A MATRIX BY A MATRIX
	Matrix matrixMulti(const Matrix& other) const;
	//--------------------------------
	//MULTIPLYING A MATRIX BY A Complex MATRIX
	cMatrix matrixMulti(const cMatrix& other) const;
	//--------------------------------
	//AUGMENTED MATRIX (MATRIX WITH VECTOR)
	Matrix augMatrix(const Vector& vector) const;
	//--------------------------------
	//AUGMENTED MATRIX (MATRIX WITH MATRIX)
	Matrix augMatrix(const Matrix& other) const;
	//--------------------------------
	//THE DETERMINANT OF A MATRIX
	float det() const;
	//--------------------------------
	//LEADING PRINCIPAL MINOR (Upper-Left 'k' by 'k' sub-matrix)
	float leadMinor(int n) const;
	//--------------------------------
	//TRANSPOSE
	Matrix transp() const;
	//--------------------------------
	//COMPARISON (Returns TRUE if current matrix is equal to Other)
	bool equalTo(const Matrix& other) const;
	//--------------------------------
	//NORM OF DIAGONAL MATRIX
	float dNorm() const;

	//--------------------------------
	//VARIABLES
	int sizeR; //the number of rows
	int sizeC; //the number of columns
	double e[MAX][MAX];
};

//====================================
//IDENTITY MATRIX
Matrix Ident(int size);

//====================================
//OUTER PRODUCT
Matrix outProduct(const Vector& u, const Vector& v);

#endif // MATRIX_H
