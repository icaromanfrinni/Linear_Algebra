#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <cmath>
#include "Matrix.h"

using namespace std;

//DEFAULT CONSTRUCTOR
Matrix::Matrix()
{
	sizeR = 0;
	sizeC = 0;
	for (int i = 0; i < MAX; i++)
		for (int j = 0; j < MAX; j++)
			e[i][j] = 0.0;
}

//OVERLOAD CONSTRUCTOR
Matrix::Matrix(int nR, int nC)
{
	sizeR = nR;
	sizeC = nC;
	for (int i = 0; i < MAX; i++)
		for (int j = 0; j < MAX; j++)
			e[i][j] = 0.0;
}


//DESTRUCTOR
Matrix::~Matrix()
{

}

//------------------------------------------------------------------------------
//PRINT
//------------------------------------------------------------------------------
void Matrix::print() const
{
	for (int i = 0; i < sizeR; i++)
	{
		for (int j = 0; j < sizeC; j++)
			cout << "\t" << e[i][j];
		cout << "\n" << endl;
	}
}

//------------------------------------------------------------------------------
//MATRIX ADDITION
//------------------------------------------------------------------------------
Matrix Matrix::matrixAdd(const Matrix &other) const
{
	Matrix R(sizeR, sizeC);
	for (int i = 0; i < sizeR; i++)
		for (int j = 0; j < sizeC; j++)
			R.e[i][j] = e[i][j] + other.e[i][j];
	return R;
}

//------------------------------------------------------------------------------
//Complex MATRIX ADDITION
//------------------------------------------------------------------------------
cMatrix Matrix::matrixAdd(const cMatrix &other) const
{
	cMatrix R(sizeR, sizeC);
	for (int i = 0; i < sizeR; i++)
		for (int j = 0; j < sizeC; j++)
			R.e[i][j] = e[i][j] + other.e[i][j];
	return R;
}

//------------------------------------------------------------------------------
//SCALAR MULTIPLICATION
//------------------------------------------------------------------------------
Matrix Matrix::scalarMulti(float factor) const
{
	Matrix R(sizeR, sizeC);
	for (int i = 0; i < sizeR; i++)
		for (int j = 0; j < sizeC; j++)
			R.e[i][j] = e[i][j] * factor;
	return R;
}

//------------------------------------------------------------------------------
//MULTIPLYING A MATRIX BY A VECTOR
//------------------------------------------------------------------------------
Vector Matrix::vectorMulti(const Vector &v) const
{
	if (sizeC != v.size)
	{
		cerr << "\n!!! INNER MATRIX DIMENSIONS MUST AGREE !!!\n" << endl;
		exit(EXIT_FAILURE);
	}
	Vector R(sizeR);
	for (int j = 0; j < sizeC; j++)
		for (int i = 0; i < sizeR; i++)
			R.e[i] += e[i][j] * v.e[j];
	return R;
}

//------------------------------------------------------------------------------
//MULTIPLYING A MATRIX BY A MATRIX
//------------------------------------------------------------------------------
Matrix Matrix::matrixMulti(const Matrix &other) const
{
	if (sizeC != other.sizeR)
	{
		cerr << "\n!!! INNER MATRIX DIMENSIONS MUST AGREE !!!\n" << endl;
		exit(EXIT_FAILURE);
	}
	Matrix R(sizeR, other.sizeC);
	for (int j = 0; j < other.sizeC; j++)
		for (int k = 0; k < sizeR; k++)
			for (int i = 0; i < sizeR; i++)
				R.e[i][j] += e[i][k] * other.e[k][j];
	return R;
}

//------------------------------------------------------------------------------
//MULTIPLYING A MATRIX BY A Complex MATRIX
//------------------------------------------------------------------------------
cMatrix Matrix::matrixMulti(const cMatrix &other) const
{
	cMatrix R(sizeR, other.sizeC);
	for (int j = 0; j < other.sizeC; j++)
		for (int k = 0; k < sizeR; k++)
			for (int i = 0; i < sizeR; i++)
				R.e[i][j] += e[i][k] * other.e[k][j];
	return R;
}

//------------------------------------------------------------------------------
//AUGMENTED MATRIX (MATRIX WITH VECTOR)
//------------------------------------------------------------------------------
Matrix Matrix::augMatrix(const Vector& vector) const
{
	Matrix aM(sizeR, sizeC + 1);
	for (int i = 0; i < sizeR; i++)
	{
		for (int j = 0; j < sizeC; j++)
			aM.e[i][j] = e[i][j];
		aM.e[i][sizeC] = vector.e[i];
	}
	return aM;
}

//------------------------------------------------------------------------------
//AUGMENTED MATRIX (MATRIX WITH MATRIX)
//------------------------------------------------------------------------------
Matrix Matrix::augMatrix(const Matrix& other) const
{
	Matrix aM(sizeR, sizeC + other.sizeC);
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
//THE DETERMINANT OF A MATRIX
//------------------------------------------------------------------------------
float Matrix::det() const
{
	if (sizeR != sizeC) //verifica se a matrix Ã© quadrada
	{
		cout << "\n!!! MATRIX MUST BE SQUARE !!!\n" << endl;
		exit(EXIT_FAILURE);
	}
	else {
		int n = sizeR; //pega a dimensÃ£o da matrix
		float det = 0.0;
		if (n == 1)
		{
			det = e[0][0];
			return det;
		}
		else if (n == 2) //se for 2x2, calcula o determinante direto
		{
			det = e[0][0] * e[1][1] - e[0][1] * e[1][0];
			return det;
		}
		else {
			for (int k = 0; k < n; k++)
			{
				Matrix subM(n - 1, n - 1); //Inicializa uma SubMatriz com DIM = n-1 da matriz corrente

				for (int i = 1; i < n; i++) //cada submatriz Ã© formada pelos elementos abaixo da primeira linha
					for (int j = 0; j < n; j++)
					{
						if (j < k) //nÃ£o pega os elementos da coluna do coeficiente
							subM.e[i - 1][j] = e[i][j];
						else if (j > k) subM.e[i - 1][j - 1] = e[i][j];
					}

				det += pow(-1, k) * e[0][k] * subM.det(); //calcula o determinante de cada submatriz e acumula
			}
			return det;
		}
	}
}

//------------------------------------------------------------------------------
//LEADING PRINCIPAL MINOR (Upper-Left 'k' by 'k' sub-matrix)
//------------------------------------------------------------------------------
float Matrix::leadMinor(int n) const
{
	Matrix subM(n, n);

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			subM.e[i][j] = e[i][j];

	return subM.det();
}

//------------------------------------------------------------------------------
//TRANSPOSE
//------------------------------------------------------------------------------
Matrix Matrix::transp() const
{
	Matrix T(sizeC, sizeR);

	for (int i = 0; i < sizeR; i++)
		for (int j = 0; j < sizeC; j++)
			T.e[j][i] = e[i][j];

	return T;
}

//------------------------------------------------------------------------------
//COMPARISON 'equal to'
//------------------------------------------------------------------------------
bool Matrix::equalTo(const Matrix &other) const
{
	if (sizeR == other.sizeR && sizeC == other.sizeC)
	{
		for (int i = 0; i < sizeR; i++)
			for (int j = 0; j < sizeC; j++)
				if (e[i][j] == other.e[i][j])
					continue;
				else return false;
		return true;
	}
	else return false;
}

//------------------------------------------------------------------------------
//IDENTITY MATRIX
//------------------------------------------------------------------------------
Matrix Ident(int size)
{
	Matrix I(size, size);
	for (int i = 0; i < size; i++)
		I.e[i][i] = 1;
	return I;
}

//------------------------------------------------------------------------------
//OUTER PRODUCT (u x v^t)
//------------------------------------------------------------------------------
Matrix outProduct(const Vector& u, const Vector& v)
{
	int m = u.size;
	int n = v.size;
	Matrix A(m, n);
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			A.e[i][j] = u.e[i] * v.e[j];
	return A;
}

//------------------------------------------------------------------------------
//THE NORM OF DIAGONAL MATRIX
//------------------------------------------------------------------------------
float Matrix::dNorm() const
{
	float N = 0.0;

	for (int i = 0; i < sizeR; i++)
		N += pow(e[i][i], 2.0);

	N = sqrt(N);
	return N;
}
