#ifndef LINEARALGEBRA_H
#define LINEARALGEBRA_H

#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <iomanip>
#include "constants.h"
#include "Matrix.h"
#include "Vector.h"
#include "cMatrix.h"
#include "cVector.h"

using namespace std;

//------------------------------------------------------------------------------
//PRINT LIKE AUGMENTED MATRIX (Matrix 'A' with Vector 'b')
//------------------------------------------------------------------------------
void printAugMatrix(const Matrix& A, const Vector& b)
{
	for (int i = 0; i < A.sizeR; i++)
	{
		for (int j = 0; j < A.sizeC; j++)
			cout << "\t" << A.e[i][j];
		cout << "\t|\t" << b.e[i] << endl;
		cout << endl;
	}
}

//------------------------------------------------------------------------------
//PRINT LIKE AUGMENTED MATRIX (Matrix 'A' with Matrix 'B')
//------------------------------------------------------------------------------
void printAugMatrix(const Matrix& A, const Matrix& B)
{
	for (int i = 0; i < A.sizeR; i++)
	{
		for (int j = 0; j < A.sizeC; j++)
			cout << "\t" << A.e[i][j];
		cout << "\t|";
		for (int j = 0; j < B.sizeC; j++)
			cout << "\t" << B.e[i][j];
		cout << "\n" << endl;
	}
}

//------------------------------------------------------------------------------
//RETORNA O ÍNDICE DO MAIOR ELEMENTO ABAIXO DA LINHA DO PIVOT
//------------------------------------------------------------------------------
int newPivot(const Matrix& A, int p)
{
	float eMax = A.e[p][p];
	int iMax = 0;
	for (int i = p + 1; i < A.sizeR; i++)
		if (A.e[i][p] > eMax)
		{
			eMax = A.e[i][p];
			iMax = i;
		}

	return iMax;
}

//------------------------------------------------------------------------------
//PARTIAL MATRIX PIVOTING (with ROW interchange)
//------------------------------------------------------------------------------
void switchRows(Matrix& M, int pivot, int nRow)
{
	float aux = 0.0;
	for (int k = 0; k < M.sizeC; k++)
	{
		aux = M.e[pivot][k];
		M.e[pivot][k] = M.e[nRow][k];
		M.e[nRow][k] = aux;
	}
}

//------------------------------------------------------------------------------
//PARTIAL VECTOR PIVOTING (with ROW interchange)
//------------------------------------------------------------------------------
void switchRows(Vector& b, int pivot, int nRow)
{
	float aux = 0.0;
	aux = b.e[pivot];
	b.e[pivot] = b.e[nRow];
	b.e[nRow] = aux;
}

//------------------------------------------------------------------------------
//PARTIAL PIVOTING (with COLUMN interchange)
//------------------------------------------------------------------------------
void switchColumns(Matrix& M, Vector& x, int pivot, int nColumn)
{
	float aux = 0.0;
	//Troca coluna da Matrix
	for (int k = 0; k < M.sizeR; k++)
	{
		aux = M.e[k][pivot];
		M.e[k][pivot] = M.e[k][nColumn];
		M.e[k][nColumn] = aux;
	}
	//Troca índice do Vetor de incógnitas
	int ind = x.flag[pivot];
	x.flag[pivot] = x.flag[nColumn];
	x.flag[nColumn] = ind;
}

//------------------------------------------------------------------------------
//FORWARD SUBSTITUTION
//------------------------------------------------------------------------------
void forSub(Matrix& L, Vector& y, const Vector& b)
{
	for (int i = 0; i < L.sizeR; i++)
	{
		if (abs(L.e[i][i]) < 0.01)
		{
			cerr << "\n\t!!! The MATRIX is Singular !!!\n" << endl;
			system("pause");
			exit(EXIT_FAILURE);
		}
		else {
			float soma = 0.0;
			for (int j = 0; j < i; j++)
				soma += L.e[i][j] * y.e[j];
			y.e[i] = (b.e[i] - soma) / L.e[i][i];
		}
	}
}

//------------------------------------------------------------------------------
//COMPLEX FORWARD SUBSTITUTION
//------------------------------------------------------------------------------
void cForSub(cMatrix& L, cVector& y, const cVector& b)
{
	for (int i = 0; i < L.sizeR; i++)
	{
		complex<double> soma = 0.0;
		for (int j = 0; j < i; j++)
			soma += L.e[i][j] * y.e[j];
		y.e[i] = (b.e[i] - soma) / L.e[i][i];
	}
}

//------------------------------------------------------------------------------
//BACK SUBSTITUTION
//------------------------------------------------------------------------------
void backSub(Matrix& U, Vector& x, Vector& c)
{
	for (int i = x.size - 1; i >= 0; i--)
	{
		if (abs(U.e[i][i]) < 0.01)
		{
			cerr << "\n\t!!! The MATRIX is Singular !!!\n" << endl;
			system("pause");
			exit(EXIT_FAILURE);
		}
		float soma = 0.0;
		for (int j = i + 1; j < U.sizeC; j++)
			soma += U.e[i][j] * x.e[x.flag[j]];
		x.e[x.flag[i]] = (c.e[i] - soma) / U.e[i][i];
	}
}

//------------------------------------------------------------------------------
//COMPLEX BACK SUBSTITUTION
//------------------------------------------------------------------------------
void cBackSub(cMatrix& U, cVector& x, cVector& c)
{
	for (int i = x.size - 1; i >= 0; i--)
	{
		complex<double> soma = 0.0;
		for (int j = i + 1; j < U.sizeC; j++)
			soma += U.e[i][j] * x.e[x.flag[j]];
		x.e[x.flag[i]] = (c.e[i] - soma) / U.e[i][i];
	}
}

//------------------------------------------------------------------------------
//GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
//------------------------------------------------------------------------------
void gaussPP(const Matrix& A, const Vector& b)
{
	cout << "\n\t*** ORIGINAL LINEAR SYSTEM ***\n" << endl;
	printAugMatrix(A, b);

	int n = A.sizeR;
	Matrix U = A;
	Vector c = b;
	Vector x(n);

	cout << "\n\t*** AUGMENTED MATRIX ***\n" << endl;
	printAugMatrix(U, c);

	//GAUSSIAN ELIMINATION
	for (int j = 0; j < n - 1; j++) //resolve cada coluna por vez
	{
		//Verifica se o 'pivot' é muito pequeno
		if (abs(U.e[j][j]) < 0.01)
		{
			float eMax = abs(U.e[j][j]); //elemento de maior valor
			int iMax = 0; //índice da linha com maior valor
			for (int i = j; i < n; i++)
				if (abs(U.e[i][j]) > eMax)
				{
					eMax = abs(U.e[i][j]);
					iMax = i;
				}

			if (iMax == 0)
			{
				cerr << "\n\t!!! The MATRIX is Singular !!!\n" << endl;
				system("pause");
				exit(EXIT_FAILURE);
			}
			else
			{
				switchRows(U, j, iMax);
				switchRows(c, j, iMax);
				cout << "\n\t!!! PARTIAL PIVOTING (switch rows) !!!\n" << endl;
				printAugMatrix(U, c);
			}
		}

		//UPPER TRIANGULAR
		for (int i = j + 1; i < U.sizeR; i++) //zera cada linha da coluna atual
		{
			float gamma = -U.e[i][j] / U.e[j][j]; //fator
			U.e[i][j] = 0.0; //zera o elemento atual
			for (int k = j + 1; k < U.sizeC; k++)
				U.e[i][k] += gamma * U.e[j][k]; //altera o resto da linha
			c.e[i] += gamma * c.e[j];
		}
	}

	//BACK SUBSTITUTION
	backSub(U, x, c);

	//OUTPUT
	cout << "\n\t*** GAUSSIAN ELIMINATION ***\n" << endl;
	printAugMatrix(U, c);

	cout << "\n\t*** SOLUTION VECTOR ***\n" << endl;
	x.print();

	cout << "\n\t======================" << endl;
	cout << "\t*** PROOF: A.x = b ***\n" << endl;
	printAugMatrix(A, x);
	Vector proof = A.vectorMulti(x);
	proof.print();
}

//------------------------------------------------------------------------------
//GAUSSIAN ELIMINATION WITH FULL PIVOTING
//------------------------------------------------------------------------------
void gaussFP(const Matrix& A, const Vector& b)
{
	cout << "\n\t*** ORIGINAL LINEAR SYSTEM ***\n" << endl;
	printAugMatrix(A, b);

	int n = A.sizeR;
	Matrix U = A;
	Vector c = b;
	Vector x(n);

	cout << "\n\t*** AUGMENTED MATRIX ***\n" << endl;
	printAugMatrix(U, c);

	//GAUSSIAN ELIMINATION
	for (int j = 0; j < n - 1; j++) //resolve cada coluna por vez
	{
		//Verifica se o 'pivot' é muito pequeno
		if (abs(U.e[j][j]) < 0.01)
		{
			float eMax = abs(U.e[j][j]); //elemento de maior valor
			int iMax = 0; //índice da linha com maior valor
			int kMax = 0; //índice da coluna com maior valor
			for (int i = j; i < n; i++)
				for (int k = j; k < n; k++)
					if (abs(U.e[i][k]) > eMax)
					{
						eMax = abs(U.e[i][k]);
						iMax = i;
						kMax = k;
					}

			if (iMax == 0 && kMax == 0)
			{
				cerr << "\n\t!!! The MATRIX is Singular !!!\n" << endl;
				system("pause");
				exit(EXIT_FAILURE);
			}
			else
			{
				switchRows(U, j, iMax);
				switchRows(c, j, iMax);
				switchColumns(U, x, j, kMax);
				cout << "\n\t!!! TOTAL PIVOTING !!!\n" << endl;
				printAugMatrix(U, c);
			}
		}

		//UPPER TRIANGULAR
		for (int i = j + 1; i < U.sizeR; i++) //zera cada linha da coluna atual
		{
			float gamma = -U.e[i][j] / U.e[j][j]; //fator
			U.e[i][j] = 0.0; //zera o elemento atual
			for (int k = j + 1; k < U.sizeC; k++)
				U.e[i][k] += gamma * U.e[j][k]; //altera o resto da linha
			c.e[i] += gamma * c.e[j];
		}
	}

	//BACK SUBSTITUTION
	backSub(U, x, c);

	//OUTPUT
	cout << "\n\t*** GAUSSIAN ELIMINATION ***\n" << endl;
	printAugMatrix(U, c);

	cout << "\n\t*** SOLUTION VECTOR ***\n" << endl;
	x.print();

	cout << "\n\t======================" << endl;
	cout << "\t*** PROOF: A.x = b ***\n" << endl;
	printAugMatrix(A, x);
	Vector proof = A.vectorMulti(x);
	proof.print();
}

//------------------------------------------------------------------------------
//GAUSS-JORDAN ELIMINATION
//------------------------------------------------------------------------------
void gaussJordan(const Matrix& A, const Vector& b)
{
	cout << "\n\t*** ORIGINAL LINEAR SYSTEM ***\n" << endl;
	printAugMatrix(A, b);

	int n = A.sizeR;
	Matrix mLeft(n, n);
	Matrix mRight(n, n);
	Vector x(n);

	cout << "\n\t*** AUGMENTED MATRIX ***\n" << endl;
	//printAugMatrix(A,Ident(n));

	//real augmented matrix
	Matrix M = A.augMatrix(Ident(n));
	M.print();

	//GAUSSIAN ELIMINATION
	for (int j = 0; j < n - 1; j++) //resolve cada coluna por vez
	{
		//Verifica se o 'pivot' é muito pequeno
		if (abs(M.e[j][j]) < 0.01)
		{
			float eMax = abs(M.e[j][j]); //elemento de maior valor
			int iMax = 0; //índice da linha com maior valor
			int kMax = 0; //índice da coluna com maior valor
			for (int i = j; i < n; i++)
				for (int k = j; k < n; k++)
					if (abs(M.e[i][k]) > eMax)
					{
						eMax = abs(M.e[i][k]);
						iMax = i;
						kMax = k;
					}

			if (iMax == 0 && kMax == 0)
			{
				cerr << "\n\t!!! The MATRIX is Singular !!!\n" << endl;
				system("pause");
				exit(EXIT_FAILURE);
			}
			else
			{
				switchRows(M, j, iMax);
				switchColumns(M, x, j, kMax);
				//cout << "\n\t!!! TOTAL PIVOTING !!!\n" << endl;
				//M.print();
			}
		}

		//UPPER TRIANGULAR
		for (int i = j + 1; i < n; i++) //zera cada linha da coluna atual
		{
			float gamma = -M.e[i][j] / M.e[j][j]; //fator
			M.e[i][j] = 0.0; //zera o elemento atual
			for (int k = j + 1; k < M.sizeC; k++)
				M.e[i][k] += gamma * M.e[j][k]; //altera o resto da linha
		}
	}

	cout << "\n\t*** UPPER TRIANGULAR ***\n" << endl;
	M.print();

	//REDUÇÃO REVERSA
	for (int j = (n - 1); j >= 0; j--)
	{
		//Verifica se o 'pivot' é muito pequeno
		if (abs(M.e[j][j]) < 0.01)
		{
			cerr << "\n\t!!! The MATRIX is Singular !!!\n" << endl;
			system("pause");
			exit(EXIT_FAILURE);
		}
		//zerar triângulo superior
		for (int i = (j - 1); i >= 0; i--)
		{
			float gamma = -M.e[i][j] / M.e[j][j];
			for (int k = 0; k < 2 * n; k++)
				M.e[i][k] += gamma * M.e[j][k];
		}
		//1 na diagonal principal
		float gamma = M.e[j][j];
		for (int k = j; k < 2 * n; k++)
		{
			if (abs(M.e[j][k]) < 0.0001)
			{
				//DEBUG: retornando zero negativo (-0.0000)
				M.e[j][k] = 0.0;
				continue;
			}
			M.e[j][k] = M.e[j][k] / gamma;
		}
	}

	//Identity Matrix
	for (int i = 0; i < M.sizeR; i++)
		for (int j = 0; j < n; j++)
			mLeft.e[i][j] = M.e[i][j];

	//Matrix Inverse
	for (int i = 0; i < M.sizeR; i++)
		for (int j = n; j < 2 * n; j++)
			mRight.e[i][j - n] = M.e[x.flag[i]][j];

	//SOLUTION VECTOR (x = A'.b)
	//Caso tenha troca de colunas, não preciso me preocupar com o array de índices
	//pq a matrix inversa já foi re-escrita na ordem correta ;P
	x = mRight.vectorMulti(b);
	/*
	Vector solVector = mRight.vectorMulti(b);
	for (int i = 0; i < x.size; i++)
		x.e[x.flag[i]] = solVector.e[i];
		*/

		//OUTPUT
	cout << "\n\t*** GAUSS-JORDAN ELIMINATION ***\n" << endl;
	printAugMatrix(mLeft, mRight);

	cout << "\n\t*** SOLUTION VECTOR ***\n" << endl;
	x.print();

	cout << "\n\t========================" << endl;
	cout << "\t*** PROOF 1: A.x = b ***\n" << endl;
	printAugMatrix(A, x);
	Vector proof1 = A.vectorMulti(x);
	proof1.print();

	cout << "\n\t=========================" << endl;
	cout << "\t*** PROOF 2: A.A' = I ***\n" << endl;
	printAugMatrix(A, mRight);
	Matrix proof2 = A.matrixMulti(mRight);
	proof2.print();
}

//------------------------------------------------------------------------------
//INVERSE
//------------------------------------------------------------------------------
Matrix inv(const Matrix& A)
{
	int n = A.sizeR;
	//Matrix mLeft(n, n);
	Matrix mRight(n, n);
	Vector x(n);

	//augmented matrix
	Matrix M = A.augMatrix(Ident(n));

	//GAUSSIAN ELIMINATION
	for (int j = 0; j < n - 1; j++) //resolve cada coluna por vez
	{
		//Verifica se o 'pivot' é muito pequeno
		if (abs(M.e[j][j]) < 0.01)
		{
			float eMax = abs(M.e[j][j]); //elemento de maior valor
			int iMax = 0; //índice da linha com maior valor
			int kMax = 0; //índice da coluna com maior valor
			for (int i = j; i < n; i++)
				for (int k = j; k < n; k++)
					if (abs(M.e[i][k]) > eMax)
					{
						eMax = abs(M.e[i][k]);
						iMax = i;
						kMax = k;
					}

			if (iMax == 0 && kMax == 0)
			{
				cerr << "\n\t!!! The MATRIX is Singular !!!\n" << endl;
				system("pause");
				exit(EXIT_FAILURE);
			}
			else
			{
				switchRows(M, j, iMax);
				switchColumns(M, x, j, kMax);
			}
		}

		//UPPER TRIANGULAR
		for (int i = j + 1; i < n; i++) //zera cada linha da coluna atual
		{
			float gamma = -M.e[i][j] / M.e[j][j]; //fator
			M.e[i][j] = 0.0; //zera o elemento atual
			for (int k = j + 1; k < M.sizeC; k++)
				M.e[i][k] += gamma * M.e[j][k]; //altera o resto da linha
		}
	}

	//REDUÇÃO REVERSA
	for (int j = (n - 1); j >= 0; j--)
	{
		//Verifica se o 'pivot' é muito pequeno
		if (abs(M.e[j][j]) < 0.01)
		{
			cerr << "\n\t!!! The MATRIX is Singular !!!\n" << endl;
			system("pause");
			exit(EXIT_FAILURE);
		}
		//zerar triângulo superior
		for (int i = (j - 1); i >= 0; i--)
		{
			float gamma = -M.e[i][j] / M.e[j][j];
			for (int k = 0; k < 2 * n; k++)
				M.e[i][k] += gamma * M.e[j][k];
		}
		//1 na diagonal principal
		float gamma = M.e[j][j];
		for (int k = j; k < 2 * n; k++)
		{
			if (abs(M.e[j][k]) < 0.0001)
			{
				//DEBUG: retornando zero negativo (-0.0000)
				M.e[j][k] = 0.0;
				continue;
			}
			M.e[j][k] = M.e[j][k] / gamma;
		}
	}

	//Matrix Inverse
	for (int i = 0; i < M.sizeR; i++)
		for (int j = n; j < 2 * n; j++)
			mRight.e[i][j - n] = M.e[x.flag[i]][j];

	return mRight;
}

//------------------------------------------------------------------------------
//INVERSE (Complex Matrix)
//------------------------------------------------------------------------------
cMatrix invC(const cMatrix& A)
{
	int n = A.sizeR;
	//Matrix mLeft(n, n);
	cMatrix mRight(n, n);
	cVector x(n);

	//augmented matrix
	cMatrix M = A.augMatrix(IdentC(n));

	//GAUSSIAN ELIMINATION
	for (int j = 0; j < n - 1; j++) //resolve cada coluna por vez
	{
		//UPPER TRIANGULAR
		for (int i = j + 1; i < n; i++) //zera cada linha da coluna atual
		{
			complex<double> gamma = -M.e[i][j] / M.e[j][j]; //fator
			M.e[i][j] = 0.0 + 0i; //zera o elemento atual
			for (int k = j + 1; k < M.sizeC; k++)
				M.e[i][k] += gamma * M.e[j][k]; //altera o resto da linha
		}
	}

	//REDUÇÃO REVERSA
	for (int j = (n - 1); j >= 0; j--)
	{
		//zerar triângulo superior
		for (int i = (j - 1); i >= 0; i--)
		{
			complex<double> gamma = -M.e[i][j] / M.e[j][j];
			for (int k = 0; k < 2 * n; k++)
				M.e[i][k] += gamma * M.e[j][k];
		}
		//1 na diagonal principal
		complex<double> gamma = M.e[j][j];
		for (int k = j; k < 2 * n; k++)
			M.e[j][k] = M.e[j][k] / gamma;
	}

	//Matrix Inverse
	for (int i = 0; i < M.sizeR; i++)
		for (int j = n; j < 2 * n; j++)
			mRight.e[i][j - n] = M.e[x.flag[i]][j];

	return mRight;
}

//------------------------------------------------------------------------------
//ORTHONORMAL BASIS VECTORS (by GRAM-SCHMIDT Process)
//------------------------------------------------------------------------------
Matrix orthonormal(const Matrix& BASE)
{
	//NÃO FAZER DETERMINANTE (É MUITO CARO)
	/*
	if (BASE.det() == 0.0)
	{
		cerr << "\n\t!!! LINEARLY DEPENDENT VECTORS (Det = 0) !!!\n" << endl;
		system("pause");
		exit(EXIT_FAILURE);
	}
	*/

	/** TRANSFORMAR A MATRIZ DE VETORES COLUNAS EM MATRIZ DE VETORES LINHA **/
	//Matrix BASE = dBASE.transp();

	//Dimension of subspace R^n
	const int n = BASE.sizeR;

	cout << "\n\t*** LINEARLY INDEPENDENT VECTORS ***\n" << endl;
	BASE.print();

	/*** STEP 1 ***/

	//Array de vetores coluna da Matrix corrente (vetores LI, i.e., BASE)
	Vector w[MAX];
	for (int j = 0; j < BASE.sizeC; j++)
		for (int i = 0; i < BASE.sizeR; i++)
			w[j].e[i] = BASE.e[i][j];

	/*** STEP 2 ***/

	//Array de vetores LI ortogonais (base ortogonal)
	Vector v[MAX];

	//First orthogonal basis element
	for (int i = 0; i < n; i++)
		v[0].e[i] = w[0].e[i];

	//Other basis vectors
	for (int j = 1; j < n; j++)
	{
		Vector soma(n);
		for (int i = 0; i < j; i++)
			soma = soma.Add(v[i].scalarMulti(w[j].Dot(v[i]) / pow(v[i].Norm(), 2.0)));
		v[j] = w[j].Sub(soma);
	}

	//Matrix BASE ORTOGONAL
	Matrix V(n, n);
	for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
			V.e[i][j] = v[j].e[i];

	cout << endl;
	cout << "\t*** ORTHOGONAL BASIS VECTORS ***" << endl;
	cout << "\t    (by GRAM-SCHMIDT Process)   " << endl;
	cout << endl;
	V.print();

	/*** STEP 3 ***/
	//Get the orthonormal basis by dividing each vector by its length

	//Array de vetores ortonormais
	Vector u[MAX];
	for (int i = 0; i < n; i++)
		u[i] = v[i].scalarMulti(1 / v[i].Norm());

	/*** RETURN ***/

	//Matrix BASE ORTONORMAL
	Matrix R(n, n);
	for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
			R.e[i][j] = u[j].e[i];

	cout << "\n\t*** ORTHONORMAL BASIS VECTORS ***\n" << endl;
	R.print();

	return R;
}

//------------------------------------------------------------------------------
//LU DECOMPOSITION
//------------------------------------------------------------------------------
void LU(const Matrix& M, const Vector& v)
{
	//Inicializar Matriz e Vetor variáveis
	Matrix A = M;
	Vector b = v;

	/* ALTO CUSTO COMPUTACIONAL */
	////Verifica se a matrix é singular, i.e., det(A) = 0
	//if (A.det() == 0.0)
	//{
	//	cerr << "\n\t!!! The MATRIX is Singular !!!\n" << endl;
	//	system("pause");
	//	exit(EXIT_FAILURE);
	//}

	cout << "\n\t*** ORIGINAL LINEAR SYSTEM ***\n" << endl;
	printAugMatrix(A, b);

	//Dimensão da matriz quadrada
	int n = b.size;
	//Inicializa Upper U = 0;
	Matrix U(n, n);
	//Inicializa Lower L = Identidade;
	Matrix L = Ident(n);
	//Inicializa vetores = 0
	Vector y(n);
	Vector x(n);

	for (int j = 0; j < n; j++)
	{
		//Verifica se o 'pivot' é muito pequeno
		if (abs(A.e[j][j]) < 0.01)
		{
			float eMax = abs(A.e[j][j]); //elemento de maior valor
			int iMax = 0; //índice da linha com maior valor
			int kMax = 0; //índice da coluna com maior valor
			for (int i = j; i < n; i++)
				for (int k = j; k < n; k++)
					if (abs(A.e[i][k]) > eMax)
					{
						eMax = abs(A.e[i][k]);
						iMax = i;
						kMax = k;
					}

			if (iMax == 0 && kMax == 0)
			{
				cerr << "\n\t!!! The MATRIX is Singular !!!\n" << endl;
				system("pause");
				exit(EXIT_FAILURE);
			}
			else
			{
				switchRows(A, j, iMax);
				switchRows(b, j, iMax);
				switchColumns(A, x, j, kMax);
				cout << "\n\t!!! TOTAL PIVOTING !!!\n" << endl;
				printAugMatrix(A, b);
			}

		}

		/*** UPPER ***/
		for (int i = 0; i <= j; i++)
		{
			float soma = 0.0;
			//for (int k = 0; k <= i - 1; k++)
			for (int k = 0; k < i; k++)
				soma += L.e[i][k] * U.e[k][j];
			U.e[i][j] = A.e[i][j] - soma;
		}
		/*** LOWER ***/
		for (int i = j + 1; i < n; i++)
		{
			float soma = 0.0;
			//for (int k = 0; k <= j - 1; k++)
			for (int k = 0; k < j; k++)
				soma += L.e[i][k] * U.e[k][j];
			L.e[i][j] = (A.e[i][j] - soma) / U.e[j][j];
		}
	}

	cout << "\n\t*** LOWER TRIANGULAR ***\n" << endl;
	L.print();
	cout << "\n\t*** UPPER TRIANGULAR ***\n" << endl;
	U.print();

	forSub(L, y, b);
	cout << "\n\t*** FORWARD SUBSTITUTION (L.y=b) ***\n" << endl;
	y.print();

	backSub(U, x, y);
	cout << "\n\t*** BACK SUBSTITUTION (U.x=y) ***\n" << endl;
	x.print();

	cout << "\n\t====================================" << endl;
	cout << "\t***        PROOF: A.x = b        ***\n" << endl;
	printAugMatrix(M, x);
	Vector proof = M.vectorMulti(x);
	proof.print();

	Matrix LU = L.matrixMulti(U);
	cout << "\n\t====================================" << endl;
	cout << "\t***        PROOF: L.U = A        ***\n" << endl;
	LU.print();
	if (LU.equalTo(M) == true)
	{
		cout << "\n\tLU DECOMPOSITION OF ORIGINAL MATRIX\n" << endl;
	}
	else {
		cout << "\n\tLU DECOMPOSITION OF PERMUTATED MATRIX\n" << endl;
	}
}

//------------------------------------------------------------------------------
//CHOLESKY DECOMPOSITION
//------------------------------------------------------------------------------
void Cholesky(const Matrix& A, const Vector& b)
{
	//Verificar se a matriz é simétrica (A == A transposto)
	if (A.equalTo(A.transp()) == false)
	{
		cerr << "\n\t!!! THE MATRIX MUST BE SYMMETRIC !!!\n" << endl;
		system("pause");
		exit(EXIT_FAILURE);
	}

	/**** NÃO FAZER A DETERMINANTE ****/
	/**** É MUITO CARO ****/
	/**** MELHOR VERIFICAR DURANTE ****/
	////Verificar se os determinantes das submatrizes são positivos (critério de Sylvester)
	//for (int k = 1; k <= A.sizeR; k++)
	//	if (A.leadMinor(k) <= 0)
	//	{
	//		cerr << "\n\t!!! THE MATRIX IS NOT POSITIVE DFINITE !!!\n" << endl;
	//		system("pause");
	//		exit(EXIT_FAILURE);
	//	}

	cout << "\n\t*** ORIGINAL LINEAR SYSTEM ***\n" << endl;
	printAugMatrix(A, b);

	int n = A.sizeR;

	//Upper triangular
	Matrix R(n, n);
	for (int i = 0; i < n; i++)
	{
		//elementos da diagonal
		float soma = 0.0;
		for (int k = 0; k < i; k++)
			soma += pow(R.e[k][i], 2.0);
		//Verifica se é positiva definida
		if (soma >= A.e[i][i])
		{
			cerr << "\n\t!!! THE MATRIX IS NOT POSITIVE DFINITE !!!\n" << endl;
			system("pause");
			exit(EXIT_FAILURE);
		}
		else R.e[i][i] = sqrt(A.e[i][i] - soma);

		//elementos acima da diagonal
		for (int j = i + 1; j < n; j++)
		{
			float soma = 0.0;
			for (int k = 0; k < i; k++)
				soma += R.e[k][i] * R.e[k][j];
			R.e[i][j] = (A.e[i][j] - soma) / R.e[i][i];
		}
	}
	//Lower triangular
	Matrix Rt = R.transp();

	cout << "\n\t*** UPPER TRIANGULAR ***\n" << endl;
	R.print();
	cout << "\n\t*** LOWER TRIANGULAR ***\n" << endl;
	Rt.print();

	Vector y(n);
	forSub(Rt, y, b);
	cout << "\n\t*** FORWARD SUBSTITUTION (Rt.y=b) ***\n" << endl;
	y.print();

	Vector x(n);
	backSub(R, x, y);
	cout << "\n\t*** BACK SUBSTITUTION (R.x=y) ***\n" << endl;
	x.print();

	cout << "\n\t====================================" << endl;
	cout << "\t***        PROOF: A.x = b        ***\n" << endl;
	printAugMatrix(A, x);
	Vector proof = A.vectorMulti(x);
	proof.print();

	Matrix RRt = Rt.matrixMulti(R);
	cout << "\n\t====================================" << endl;
	cout << "\t***       PROOF: Rt.R = A        ***\n" << endl;
	RRt.print();
}

//------------------------------------------------------------------------------
//REDUCED ROW ECHELON FORM
//------------------------------------------------------------------------------
Matrix RREF(const Matrix& M)
{
	Matrix A = M;
	cout << "\n\t*** MATRIX ***\n" << endl;
	A.print();

	//FORWARD PHASE
	for (int i = 0; i < A.sizeR; i++)
	{
		double pivot = 0.0;
		int col = 0;
		for (int j = 0; j < A.sizeC; j++)
			if (A.e[i][j] != 0.0)
			{
				pivot = A.e[i][j];
				col = j;
				break;
			}

		if (pivot != 0.0)
		{
			for (int j = col; j < A.sizeC; j++)
			{
				if (abs(A.e[i][j]) < 0.0001)
				{
					//DEBUG: retornando zero negativo (-0.00)
					A.e[i][j] = 0.0;
					continue;
				}
				A.e[i][j] = A.e[i][j] / pivot;
			}
			for (int k = i + 1; k < A.sizeR; k++)
			{
				double gamma = -A.e[k][col];

				for (int j = 0; j < A.sizeC; j++)
				{
					//cout << "\tgama = " << gamma << endl;
					//cout << "\tA.e[" << i + 1 << "][" << j + 1 << "] = " << A.e[i][j] << endl;
					A.e[k][j] += gamma * A.e[i][j];
					//cout << "\tA.e[" << k + 1 << "][" << j + 1 << "] = " << A.e[k][j] << endl;
				}
			}
		}
	}

	cout << "\n\t*** FORWARD PHASE ***\n" << endl;
	A.print();

	//BACKWARD PHASE
	for (int i = A.sizeR - 1; i >= 0; i--)
	{
		int col = 0;
		for (int j = 0; j < A.sizeC; j++)
			if (A.e[i][j] == 1.0)
			{
				col = j;
				break;
			}
		for (int k = 0; k < i; k++)
		{
			float gamma = -A.e[k][col];
			for (int j = 0; j < A.sizeC; j++)
				A.e[k][j] += gamma * A.e[i][j];
		}
	}
	cout << "\n\t*** BACKWARD PHASE ***\n" << endl;
	A.print();

	//LINEARLY INDEPENDENT VECTORS (Column)
	//Inicializa como uma matriz Identidade
	//para preencher com base canônica as colunas dos vetores dependentes
	Matrix LI = Ident(M.sizeR);

	//ROW INTERCHANGE (REORDER)
	int row = 0; //primeira linha de busca do '1' para cada coluna
	for (int j = 0; j < A.sizeC; j++)
	{
		for (int i = row; i < A.sizeR; i++)
			if (A.e[i][j] == 1.0)   //se achar, troca a linha corrente do '1' pela linha do início da busca (row)
			{                       //segue para a próxima coluna, começando na próxima linha abaixo do 'row'
				for (int k = 0; k < A.sizeC; k++)
				{
					float aux = A.e[i][k];
					A.e[i][k] = A.e[row][k];
					A.e[row][k] = aux;
				}

				//Matriz de Vetores Coluna LI (recebe coluna LI da Matriz original)
				for (int k = 0; k < M.sizeR; k++)
				{
					LI.e[k][row] = M.e[k][j];
				}

				row++;
				break;
			}
	}

	//OUTPUT
	cout << "\n\t*** REDUCED ROW ECHELON FORM ***\n" << endl;
	A.print();

	//orthonormal(LI);

	return A;
}

//------------------------------------------------------------------------------
//Complex REDUCED ROW ECHELON FORM
//------------------------------------------------------------------------------
cMatrix cRREF(const cMatrix& M)
{
	cMatrix A = M;
	//cout << "\n\t*** MATRIX ***\n" << endl;
	//A.print();

	//FORWARD PHASE
	for (int i = 0; i < A.sizeR; i++)
	{
		complex<double> pivot = 0.0;
		int col = 0;
		for (int j = 0; j < A.sizeC; j++)
			if (A.e[i][j] != 0.0)
			{
				pivot = A.e[i][j];
				col = j;
				break;
			}

		if (pivot != 0.0)
		{
			for (int j = col; j < A.sizeC; j++)
			{
				if (abs(A.e[i][j]) < 0.001)
				{
					//DEBUG: retornando zero negativo (-0.00)
					A.e[i][j] = 0.0;
					continue;
				}
				A.e[i][j] = A.e[i][j] / pivot;
			}
			for (int k = i + 1; k < A.sizeR; k++)
			{
				complex<double> gamma = -A.e[k][col];

				for (int j = 0; j < A.sizeC; j++)
				{
					//cout << "\tgama = " << gamma << endl;
					//cout << "\tA.e[" << i + 1 << "][" << j + 1 << "] = " << A.e[i][j] << endl;
					A.e[k][j] += gamma * A.e[i][j];
					//cout << "\tA.e[" << k + 1 << "][" << j + 1 << "] = " << A.e[k][j] << endl;
				}
			}
		}
	}

	//cout << "\n\t*** FORWARD PHASE ***\n" << endl;
	//A.print();

	//BACKWARD PHASE
	//for (int i = A.sizeR - 1; i >= 0; i--)
	//{
	//	int col = 0;
	//	for (int j = 0; j < A.sizeC; j++)
	//		if (A.e[i][j] == 1.0)
	//		{
	//			col = j;
	//			break;
	//		}
	//	for (int k = 0; k < i; k++)
	//	{
	//		complex<double> gamma = -A.e[k][col];
	//		for (int j = 0; j < A.sizeC; j++)
	//			A.e[k][j] += gamma * A.e[i][j];
	//	}
	//}
	
	//cout << "\n\t*** BACKWARD PHASE ***\n" << endl;
	//A.print();

	//OUTPUT
	cout << "\n\t*** REDUCED ROW ECHELON FORM ***\n" << endl;
	A.print();

	return A;
}

//------------------------------------------------------------------------------
//QR Decomposition (by Gram-Schmidt process)
//------------------------------------------------------------------------------
void QRbyGS(const Matrix& A, const Vector& b)
{
	//ORTHONORMAL BASIS
	Matrix Q = orthonormal(A); //*********** Matriz de vetores independentes
	//Q.print();

	//UPPER TRIANGULAR MATRIX
	Matrix Qt = Q.transp();
	Matrix R = Qt.matrixMulti(A);

	int n = A.sizeR;
	//DEBUG: retornando zero negativo (-0.0000)
	for (int j = 0; j < n - 1; j++)
		for (int i = j + 1; i < n; i++)
			if (abs(R.e[i][j]) < 0.0001) R.e[i][j] = 0.0;

	cout << "\n\t*** UPPER TRIANGULAR MATRIX ***" << endl;
	cout << "\t***  from QR Decomposition  ***" << endl;
	cout << "\t***(by Gram-Schmidt process)***\n" << endl;
	R.print();

	//Vector y = Qt * b
	Vector y = Qt.vectorMulti(b);

	//Solution vector x
	Vector x(b.size);
	backSub(R, x, y);
	cout << "\n\t*** BACK SUBSTITUTION (R.x=y) ***\n" << endl;
	x.print();

	cout << "\n\t======================================" << endl;
	cout << "\t***         PROOF: A.x = b         ***\n" << endl;
	printAugMatrix(A, x);
	Vector proof = A.vectorMulti(x);
	proof.print();

	Matrix QR = Q.matrixMulti(R);
	cout << "\n\t======================================" << endl;
	cout << "\t***         PROOF: Q.R = A         ***\n" << endl;
	QR.print();
}

//------------------------------------------------------------------------------
//QR Decomposition (by Reflection matrix or Householder)
//------------------------------------------------------------------------------
void QRbyRef(const Matrix& A, const Vector& b)
{
	// *************************** MATRIZ QUADRADA
	int n = A.sizeR;
	Matrix Qt = Ident(n);
	Matrix R = A;

	for (int j = 0; j < n - 1; j++)
	{
		Vector v(n);

		for (int i = j; i < n; i++)
			v.e[i] = R.e[i][j]; //Vetor da coluna 'j' da matriz 'A' com elementos acima zerados

		Matrix H(n, n);
		if (v.Norm() == abs(v.e[j]))
			H = Ident(n);
		else {
			Vector w(n);
			w.e[j] = v.Norm(); //Norma euclidiana de 'v' vezes base canônica de 'j'
			Vector u = v.Sub(w); //vetor perpendicular ao eixo de reflexão
			u.divScalar(u.Norm());
			Matrix U = outProduct(u, u);
			H = Ident(n).matrixAdd(U.scalarMulti(-2.0));
		}

		Qt = H.matrixMulti(Qt);
		R = Qt.matrixMulti(A);
	}

	//DEBUG: retornando zero negativo (-0.0000)
	for (int j = 0; j < n - 1; j++)
		for (int i = j + 1; i < n; i++)
			if (abs(R.e[i][j]) < 0.0001) R.e[i][j] = 0.0;

	cout << "\n\t*** UPPER TRIANGULAR MATRIX ***" << endl;
	cout << "\t***  from QR Decomposition  ***" << endl;
	cout << "\t***  (by Reflection matrix) ***\n" << endl;
	R.print();

	//Vector y = Qt * b
	Vector y = Qt.vectorMulti(b);

	//Solution vector x
	Vector x(b.size);
	backSub(R, x, y);
	cout << "\n\t*** BACK SUBSTITUTION (R.x=y) ***\n" << endl;
	x.print();

	cout << "\n\t======================================" << endl;
	cout << "\t***         PROOF: A.x = b         ***\n" << endl;
	printAugMatrix(A, x);
	Vector proof = A.vectorMulti(x);
	proof.print();

	Matrix Q = Qt.transp();
	Matrix QR = Q.matrixMulti(R);
	cout << "\n\t======================================" << endl;
	cout << "\t***         PROOF: Q.R = A         ***\n" << endl;
	QR.print();
}

//------------------------------------------------------------------------------
//HOUSEHOLDER Matrix (Eigen-value and Eigen-vector)
//------------------------------------------------------------------------------
Matrix Hj(const Matrix& A, const int& j)
{
	int n = A.sizeR;
	Vector p(n);
	for (int i = j + 1; i < n; i++) //(j + 1) pq é abaixo da diagonal
		p.e[i] = A.e[i][j]; //Vetor da coluna 'j' da matriz 'A' com elementos acima zerados

	Matrix H(n, n);
	if (p.Norm() == abs(p.e[j + 1]))
		H = Ident(n);
	else {
		Vector pL(n);
		pL.e[j + 1] = p.Norm(); //Norma euclidiana de 'v' vezes base canônica de 'j'
		Vector u = p.Sub(pL); //vetor perpendicular ao eixo de reflexão
		u.divScalar(u.Norm());
		Matrix U = outProduct(u, u);
		H = Ident(n).matrixAdd(U.scalarMulti(-2.0));
	}

	return H;
}

//------------------------------------------------------------------------------
//QR Decomposition (by Rotation matrix or Givens)
//------------------------------------------------------------------------------
void QRbyRot(const Matrix& A, const Vector& b)
{
	// ******************************* MATRIX RETANGULAR
	int m = A.sizeR;
	int n = A.sizeC;
	Matrix Qt = Ident(m);
	Matrix R = A;

	for (int j = 0; j < n - 1; j++)
		for (int i = j + 1; i < m; i++)
		{
			float teta = 0.0;
			if (abs(R.e[j][j]) < 0.0001)
				teta = PI / 2;
			else teta = atan(R.e[i][j] / R.e[j][j]);

			//Matriz de rotação
			Matrix G = Ident(m);
			G.e[i][i] = cos(teta);
			G.e[j][j] = cos(teta);
			G.e[j][i] = sin(teta);
			G.e[i][j] = -sin(teta);

			Qt = G.matrixMulti(Qt);
			R = Qt.matrixMulti(A);
		}

	//DEBUG: retornando zero negativo (-0.0000)
	for (int j = 0; j < n; j++)
		for (int i = j + 1; i < m; i++)
			if (abs(R.e[i][j]) < 0.0001) R.e[i][j] = 0.0;

	cout << "\n\t*** UPPER TRIANGULAR MATRIX ***" << endl;
	cout << "\t***  from QR Decomposition  ***" << endl;
	cout << "\t***   (by Rotation matrix)  ***\n" << endl;
	R.print();

	//Vector y = Qt * b
	Vector y = Qt.vectorMulti(b);

	//Solution vector x
	Vector x(b.size);
	backSub(R, x, y);
	cout << "\n\t*** BACK SUBSTITUTION (R.x=y) ***\n" << endl;
	x.print();

	cout << "\n\t======================================" << endl;
	cout << "\t***         PROOF: A.x = b         ***\n" << endl;
	printAugMatrix(A, x);
	Vector proof = A.vectorMulti(x);
	proof.print();

	Matrix Q = Qt.transp();
	Matrix QR = Q.matrixMulti(R);
	cout << "\n\t======================================" << endl;
	cout << "\t***         PROOF: Q.R = A         ***\n" << endl;
	QR.print();
}

//------------------------------------------------------------------------------
//	EIGEN-VALUES and EIGEN-VECTORS (by Power Method)
//
//	BIBLIOGRAPHY:
//	Golub - Matrix Computations (2013) pg366
//	Watkins - Fundamentals of Matrix Computations (2010) pg314
//------------------------------------------------------------------------------
void powerMethod(const Matrix& A, Vector& x)
{
	cout << "\n\t*** ORIGINAL MATRIX ***\n" << endl;
	A.print();

	//Precisão
	float eps;
	cout << "\n\tEnter the accuracy desired: ";
	cin >> eps;

	//Auto-Valor inicial
	float lambda = x.e[0];
	float lambda0;

	//Auto-Vetor
	Vector q;

	//Contador
	int nCount = 0;

	do
	{
		nCount++;
		lambda0 = lambda;

		x.divScalar(x.Norm());
		q = x;
		x = A.vectorMulti(q);
		lambda = q.Dot(x);

	} while (fabs((lambda - lambda0) / lambda) > eps);

	cout << "\n\t*** THE POWER METHOD ***" << endl;

	cout << "\n\tNumber of interations: " << nCount << endl;
	cout << "\n\tEigen-value = " << lambda << endl;
	cout << "\n\tEigen-vector:\n" << endl;
	q.print();

	cout << "\n\t=================================" << endl;
	cout << "\t***   PROOF: A.v = lambda.v   ***\n" << endl;

	Vector y;
	cout << "\tMultiplying matrix 'A' by vector 'v':\n" << endl;
	y = A.vectorMulti(q);
	y.print();
	cout << "\tMultiplying scalar 'lambda' by vector 'v':\n" << endl;
	y = q.scalarMulti(lambda);
	y.print();
}

//------------------------------------------------------------------------------
//	COMPLEX EIGEN-VALUES and EIGEN-VECTORS (by Power Method)
//------------------------------------------------------------------------------
void cPowerMethod(const cMatrix& A)
{
	//Precisão
	double eps;
	cout << "\n\tEnter the accuracy desired: ";
	cin >> eps;

	//Chute inicial
	cVector x(A.sizeR);
	x.e[0] = 1.0;
	//Auto-Valor inicial
	complex<double> lambda = x.e[0];
	complex<double> lambda0;

	//Auto-Vetor
	cVector q;

	//Contador
	int nCount = 0;

	do
	{
		nCount++;
		lambda0 = lambda;

		x.divScalar(x.Norm());
		q = x;
		x = A.vectorMulti(q);
		lambda = q.Dot(x);

	} while (abs((lambda - lambda0) / lambda) > eps);

	cout << "\n\t*** THE POWER METHOD ***" << endl;

	cout << "\n\tNumber of interations: " << nCount << endl;
	cout << "\n\tEigen-value = " << lambda << endl;
	cout << "\n\tEigen-vector:\n" << endl;
	q.print();

	cout << "\n\t=================================" << endl;
	cout << "\t***   PROOF: A.v = lambda.v   ***\n" << endl;

	cVector y;
	cout << "\tMultiplying matrix 'A' by vector 'v':\n" << endl;
	y = A.vectorMulti(q);
	y.print();
	cout << "\tMultiplying scalar 'lambda' by vector 'v':\n" << endl;
	y = q.scalarMulti(lambda);
	y.print();
}

//------------------------------------------------------------------------------
//	EIGEN-VALUES and EIGEN-VECTORS (by Inverse Power Method)
//------------------------------------------------------------------------------
void invPowerMethod(const Matrix& A, Vector& x)
{
	cout << "\n\t*** ORIGINAL MATRIX ***\n" << endl;
	A.print();

	//Precisão
	float eps;
	cout << "\n\tEnter the accuracy desired: ";
	cin >> eps;

	//Auto-Valor inicial
	float lambda = x.e[0];
	float lambda0;

	//Auto-Vetor
	Vector q;

	//Contador
	int nCount = 0;

	//Matriz inversa: O(n3)
	Matrix Ai = inv(A);

	do
	{
		nCount++;
		lambda0 = lambda;

		x.divScalar(x.Norm());
		q = x;
		x = Ai.vectorMulti(q);
		lambda = q.Dot(x);

	} while (fabs((lambda - lambda0) / lambda) > eps);

	//Autovalor inverso
	lambda = 1 / lambda;

	cout << "\n\t*** THE INVERSE POWER METHOD ***" << endl;

	cout << "\n\tNumber of interations: " << nCount << endl;
	cout << "\n\tEigen-value = " << lambda << endl;
	cout << "\n\tEigen-vector:\n" << endl;
	q.print();

	cout << "\n\t=================================" << endl;
	cout << "\t***   PROOF: A.v = lambda.v   ***\n" << endl;

	Vector y;
	cout << "\tMultiplying matrix 'A' by vector 'v':\n" << endl;
	y = A.vectorMulti(q);
	y.print();
	cout << "\tMultiplying scalar 'lambda' by vector 'v':\n" << endl;
	y = q.scalarMulti(lambda);
	y.print();
}

//------------------------------------------------------------------------------
//	EIGEN-VALUES and EIGEN-VECTORS (by Inverse Power Method with Displacement)
//------------------------------------------------------------------------------
void invPowMethodDispl(const Matrix& A, Vector& x)
{
	cout << "\n\t*** ORIGINAL MATRIX ***\n" << endl;
	A.print();

	//Deslocamento
	float alpha;
	cout << "\n\tEnter the displacement: ";
	cin >> alpha;

	//Matriz de deslocamento
	Matrix D = Ident(A.sizeR);
	D = D.scalarMulti(-alpha);

	//Matriz modificada
	Matrix Am = A.matrixAdd(D);

	/* MÉTODO DA POTÊNCIA INVERSA */
	//Precisão
	float eps;
	cout << "\n\tEnter the accuracy desired: ";
	cin >> eps;

	//Auto-Valor inicial
	float lambda = x.e[0];
	float lambda0;

	//Auto-Vetor
	Vector q;

	//Contador
	int nCount = 0;

	//Matriz inversa: O(n3)
	Matrix Ai = inv(Am);

	/********* SE MATRIZ 'Am' é singular, então o deslocamento inserido é um auto-valor ************/

	do
	{
		nCount++;
		lambda0 = lambda;

		x.divScalar(x.Norm());
		q = x;
		x = Ai.vectorMulti(q);
		lambda = q.Dot(x);

	} while (fabs((lambda - lambda0) / lambda) > eps);

	//Autovalor deslocado inverso
	lambda = 1 / lambda;
	//Autovalor mais próximo
	lambda = lambda + alpha;

	cout << "\n\t*** THE INVERSE POWER METHOD ***" << endl;
	cout << "\t***     with Displacement    ***" << endl;

	cout << "\n\tNumber of interations: " << nCount << endl;
	cout << "\n\tEigen-value = " << lambda << endl;
	cout << "\n\tEigen-vector:\n" << endl;
	q.print();

	cout << "\n\t=================================" << endl;
	cout << "\t***   PROOF: A.v = lambda.v   ***\n" << endl;

	Vector y;
	cout << "\tMultiplying matrix 'A' by vector 'v':\n" << endl;
	y = A.vectorMulti(q);
	y.print();
	cout << "\tMultiplying scalar 'lambda' by vector 'v':\n" << endl;
	y = q.scalarMulti(lambda);
	y.print();
}

//------------------------------------------------------------------------------
//	EIGEN-VECTORS (of Upper Triangular Matrix - Schur form)
//------------------------------------------------------------------------------
Matrix EigenVectors(const Matrix& T)
{
	int n = T.sizeR;
	Matrix V = Ident(n);

	for (int j = 1; j < n; j++)
		for (int i = j - 1; i >= 0; i--)
		{
			float soma = 0.0;

			for (int k = i + 1; k <= j; k++)
				soma += T.e[i][k] * V.e[k][j];

			V.e[i][j] = -soma / (T.e[i][i] - T.e[j][j]);
		}

	return V;
}

//------------------------------------------------------------------------------
//	COMPLEX EIGEN-VALUES (by 2x2 Block Matrix from Hessenberg Form)
//------------------------------------------------------------------------------
cMatrix ComplexEigenValeus(const Matrix& M)
{
	int n = M.sizeR;
	cMatrix cM(n, n);
	double a, b, c, d;
	double A, B, C;
	complex<double> delta, lambda1, lambda2;

	for (int i = 0; i < n; i++)
	{
		//se matriz dentada tiver dimensão ímpar
		if (i == n - 1)
		{
			cM.e[i][i] = M.e[i][i];
			break;
		}

		//raízes do polinômio característico
		a = M.e[i][i]; b = M.e[i][i + 1]; c = M.e[i + 1][i]; d = M.e[i + 1][i + 1];
		A = 1.0; B = -(a + d); C = a * d - c * b;
		delta = pow(B, 2.0) - 4.0 * A*C;
		lambda1 = (-B + sqrt(delta)) / (2.0 * A);
		lambda2 = (-B - sqrt(delta)) / (2.0 * A);

		//matriz diagonal com as raízes (autovalores)
		cM.e[i][i] = lambda1;
		cM.e[i + 1][i + 1] = lambda2;

		i++;
	}

	return cM;
}

//------------------------------------------------------------------------------
//COMPLEX LU DECOMPOSITION (without pitoving)
//------------------------------------------------------------------------------
cVector cLU(const cMatrix& M, const cVector& v)
{
	//Inicializar Matriz e Vetor variáveis
	cMatrix A = M;
	cVector b = v;

	//Dimensão da matriz quadrada
	int n = b.size;
	//Inicializa Upper U = 0;
	cMatrix U(n, n);
	//Inicializa Lower L = Identidade;
	cMatrix L = IdentC(n);
	//Inicializa vetores = 0
	cVector y(n);
	cVector x(n);

	for (int j = 0; j < n; j++)
	{
		/*** UPPER ***/
		for (int i = 0; i <= j; i++)
		{
			complex<double> soma = 0.0;
			//for (int k = 0; k <= i - 1; k++)
			for (int k = 0; k < i; k++)
				soma += L.e[i][k] * U.e[k][j];
			U.e[i][j] = A.e[i][j] - soma;
		}
		/*** LOWER ***/
		for (int i = j + 1; i < n; i++)
		{
			complex<double> soma = 0.0;
			//for (int k = 0; k <= j - 1; k++)
			for (int k = 0; k < j; k++)
				soma += L.e[i][k] * U.e[k][j];
			L.e[i][j] = (A.e[i][j] - soma) / U.e[j][j];
		}
	}

	cout << "\n\t*** LOWER TRIANGULAR ***\n" << endl;
	L.print();
	cout << "\n\t*** UPPER TRIANGULAR ***\n" << endl;
	U.print();

	cForSub(L, y, b);
	cout << "\n\t*** FORWARD SUBSTITUTION (L.y=b) ***\n" << endl;
	y.print();

	cBackSub(U, x, y);
	cout << "\n\t*** BACK SUBSTITUTION (U.x=y) ***\n" << endl;
	x.print();

	/*
	cout << "\n\t====================================" << endl;
	cout << "\t***        PROOF: A.x = b        ***\n" << endl;
	cVector proof = M.vectorMulti(x);
	proof.print();

	cMatrix LU = L.matrixMulti(U);
	cout << "\n\t====================================" << endl;
	cout << "\t***        PROOF: L.U = A        ***\n" << endl;
	LU.print();
	*/

	return x;
}

//------------------------------------------------------------------------------
//	COMPLEX EIGEN-VECTORS (from Hessenberg Form)
//  Upper Block-Triangle Matrix
//  http://www.wright.edu/~chaocheng.huang/lecture/mth255/mth255lect1.pdf
//------------------------------------------------------------------------------
cMatrix ComplexEigenVectors(const Matrix& M, const cMatrix& cD)
{
	//PASSOS:
	//1. Matriz blocada real (M)
	//2. Para cada autovalor (lambda)
	//  a) C=M-lambda.I
	//  b) Forma escalonada reduzida (R) de C
	//  c) Autovetor (x) do lambda escolhido
	//  d) Resolver o problema: R.x=0

	int n = M.sizeR;
	cMatrix cV(n, n);

	//PARA CADA AUTOVALOR
	for (int k = 0; k < n; k++)
	{
		complex<double> lambda = cD.e[k][k];
		//C=M-lambda.I
		cMatrix C = M.matrixAdd(IdentC(n).scalarMulti(-lambda));
		//Forma escalonada reduzida
		cMatrix R(n, n);
		R = cRREF(C);

		//DEBUG
		/*if (k == 0)
		{
			R.e[0][0] = 1.0 + 0i;	R.e[0][1] = -1.21086 - 0.96103i;	R.e[0][2] = 0.0 + 0i;
			R.e[1][0] = 0.0 + 0i;	R.e[1][1] = 0.0 + 0i;				R.e[1][2] = 1.0 + 0i;
			R.e[2][0] = 0.0 + 0i;	R.e[2][1] = 0.0 + 0i;				R.e[2][2] = 0.0 + 0i;
		}
		if (k == 1)
		{
			R.e[0][0] = 1.0 + 0i;	R.e[0][1] = -1.21086 + 0.96103i;	R.e[0][2] = 0.0 + 0i;
			R.e[1][0] = 0.0 + 0i;	R.e[1][1] = 0.0 + 0i;				R.e[1][2] = 1.0 + 0i;
			R.e[2][0] = 0.0 + 0i;	R.e[2][1] = 0.0 + 0i;				R.e[2][2] = 0.0 + 0i;
		}
		if (k == 2)
		{
			R.e[0][0] = 1.0 + 0i;	R.e[0][1] = 0.0 + 0i;	R.e[0][2] = -0.28551 + 0i;
			R.e[1][0] = 0.0 + 0i;	R.e[1][1] = 1.0 + 0i;	R.e[1][2] = 0.95740 + 0i;
			R.e[2][0] = 0.0 + 0i;	R.e[2][1] = 0.0 + 0i;	R.e[2][2] = 0.0 + 0i;
		}*/

		//Autovetor do lambda escolhido
		//Se pivot é zero (Uii = 0):
		//  x[i+1] = 0
		//  x[i] = 1
		//  continua...
		//Senão 'faz a conta'
		//  soma = x2U12 + x3U13 + x4U14
		//  x1 = -soma
		cVector q(n);
		for (int i = n - 1; i > -1; i--)
		{
			if (R.e[i][i].real() < 0.001)
			{
				q.e[i + 1] = 0.0;
				q.e[i] = 1.0;
				continue;
			}
			else
			{
				complex<double> soma = 0.0;
				for (int j = i + 1; j < n; j++)
					soma += q.e[j] * R.e[i][j];
				q.e[i] = -soma;
			}
		}
		q.print();
		q.divScalar(q.Norm()); //orthonomalize
		q.print();

		//Armazena cada autovetor na respectiva coluna da matriz de autovetores
		for (int i = 0; i < n; i++)
			cV.e[i][k] = q.e[i];
	}

	return cV;
}

//------------------------------------------------------------------------------
//	EIGEN-VALUES and EIGEN-VECTORS (by Householder Method)
//------------------------------------------------------------------------------
void HouseholderMethod(const Matrix& A)
{
	int n = A.sizeR;
	Matrix H = Ident(n);

	Matrix T = A;
	for (int k = 0; k < n; k++)
	{
		Matrix Hk = Hj(T, k);
		T = Hk.matrixMulti(T.matrixMulti(Hk)); //Matriz transformada
		H = H.matrixMulti(Hk); //Acumulação da transformação
	}

	//DEBUG: retornando zero negativo (-0.0000)
	/*
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			if (abs(T.e[i][j]) < 0.0001) T.e[i][j] = 0.0;
	*/

	cout << "\n\t*** ORIGINAL MATRIX ***\n" << endl;
	A.print();

	/* DECOMPOSIÇÃO QR */

	//Inicializa a Norma da Diagonal
	float normNew = T.dNorm();
	float normOld;
	//Precisão
	float eps;
	cout << "\n\tEnter the accuracy desired: ";
	cin >> eps;
	//Inicializa um contador
	int nCount = 0;

	Matrix D = T; //Inicializa como Tridiagonal para preservar a matriz T
	Matrix Q = Ident(n); //Inicializa como Identidade para acumular
	do
	{
		nCount++;
		normOld = normNew;
		//Decomposição QR (por Gram-Schmidt)
		Q = Q.matrixMulti(orthonormal(D)); //Acumulação de transformações
		Matrix Qt = Q.transp();
		D = Qt.matrixMulti(T.matrixMulti(Q));

		normNew = D.dNorm();
	} while (fabs(normNew - normOld) > eps);

	//Inicializa Matriz de AutoVetores como Identidade
	Matrix V = Ident(n);

	//Verifica se a Matrix é simétrica
	if (A.equalTo(A.transp()) == true)
	{
		cout << "\n\t*** HOUSEHOLDER MATRIX ***\n" << endl;
		H.print();
		cout << "\n\t*** TRIDIAGONAL FORM ***\n" << endl;
		T.print();
		cout << "\n\t*** DIAGONAL MATRIX ***" << endl;
		cout << "\t***  (Eigenvaleus)  ***\n" << endl;
		D.print(); //já são os auto-valores
		V = H.matrixMulti(Q); //Autovetores (resultado da acumulação das transformações)
		cout << "\n\t*** EIGEN-VECTORS ***\n" << endl;
		V.print();

		//PROOF
		cout << "\n\t=================================" << endl;
		cout << "\t***   PROOF: A = V.D.V-1   ***\n" << endl;
		Matrix P = V.matrixMulti(D.matrixMulti(inv(V)));
		P.print();
	}
	else
	{
		cout << "\n\t*** UPPER HESSENBERG FORM ***\n" << endl;
		T.print();
		cout << "\n\t*** UPPER TRIANGULAR MATRIX ***\n" << endl;
		D.print();

		//Verifica se todos os elementos abaixo da diagonal são zeros (ou muito próximos)
		for (int i = 0; i < n - 1; i++)
			if (fabs(D.e[i + 1][i]) > 0.1)
			{
				cMatrix cD(n, n);
				cD = ComplexEigenValeus(D);
				cout << "\n\t*** EIGEN-VALEUS ***" << endl;
				cout << "\t***  (Diagonal)  ***\n" << endl;
				cD.print();

				cMatrix cV(n, n);
				cV = ComplexEigenVectors(D, cD);
				cout << "\n\t*** EIGEN-VECTORS ***\n" << endl;
				cV = Q.matrixMulti(cV);
				cV = H.matrixMulti(cV);
				cV.print();

				/* FAZER PROVA SEM Matriz Inversa -> A.V = V.D */

				//PROOF
				cout << "\n\t=================================" << endl;
				cout <<   "\t***    PROOF: A = V.D.V-1    ***\n" << endl;
				cMatrix P = cV.matrixMulti(cD.matrixMulti(invC(cV)));
				P.print();

				break;
			}
			else
			{
				/* APENAS PARA AUTOVALORES DIFERENTES */
				V = EigenVectors(D);

				//Multiplicar à esquerda a acumulação das transformações
				//para obter os autovetores da matriz original
				V = Q.matrixMulti(V); //Transformações de similaridade para obter a forma Triangular Superior
				V = H.matrixMulti(V); //Método de Householder para obter a forma de Hessenberg

				cout << "\n\t*** EIGEN-VALEUS ***" << endl;
				cout << "\t***  (Diagonal)  ***\n" << endl;
				for (int i = 0; i < n - 1; i++)
					for (int j = i + 1; j < n; j++)
						D.e[i][j] = 0;
				D.print();

				cout << "\n\t*** EIGEN-VECTORS ***\n" << endl;
				V.print();

				//PROOF
				cout << "\n\t=================================" << endl;
				cout << "\t***   PROOF: A = V.D.V-1   ***\n" << endl;
				Matrix P = V.matrixMulti(D.matrixMulti(inv(V)));
				P.print();

				break;
			}
	}

	cout << "\n\tNumber of interations: " << nCount << endl;
}


#endif // LINEARALGEBRA_H
#pragma once
