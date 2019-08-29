#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include "LinearAlgebra.h"

using namespace std;

//------------------------------------------------------------------------------
//FUNCTION DECLARATIONS
//------------------------------------------------------------------------------

//SCREEN
void clearScreen();
void pauseScreen();

//MAIN MENU
void mainMenu();
//Read Input File
bool readFile(Matrix& A, Vector& b);
//1. Multiplying a Matrix by a Vector
void multiMatrixVector();
//2. Multiplying a Matrix by a Matrix
void multiMatrixMatrix();
//3. Gaussian Elimination (Partial Pivoting)
void menuGaussPP();
//4. Gaussian Elimination (Full Pivoting)
void menuGaussFP();
//5. Gauss-Jordan Elimination
void menuGaussJordan();
//6. LU Decomposition
void menuLU();
//7. Cholesky Decomposition
void menuCholesky();
//8. Reduced Row Echelon Form
void menuRREF();
//9. The Gram-Schmidt Process
void menuOrth();
//10. QR Decomposition (by Gram-Schmidt process)
void menuQRbyGS();
//11. QR Decomposition (by Reflection matrix)
void menuQRbyRef();
//12. QR Decomposition (by Rotation matrix)
void menuQRbyRot();
//13. Eigen-value and Eigen-vector (by Power Method)
void menuPowerMethod();
//14. Eigen-value and Eigen-vector (by Inverse Power Method)
void menuInvPowerMethod();
//15. Eigen-value and Eigen-vector (by Inverse Power Method with displacement)
void menuInvPowMethodDispl();
//16. Eigen-value and Eigen-vector (by Householder Method)
void menuHouseholderMethod();

//------------------------------------------------------------------------------
//MAIN FUNCTION
//------------------------------------------------------------------------------
int main()
{
	cout << setprecision(2) << fixed;

	while (true)
		mainMenu();

	return 0;
}

//------------------------------------------------------------------------------
//FUNCTION DEFINITIONS
//------------------------------------------------------------------------------

//SCREEN
//Clear
void clearScreen()
{
#ifdef __unix__
	system("clear");
#elif defined _WIN32
	system("cls");
#endif
}
//Pause
void pauseScreen()
{
#ifdef __unix__
	system("read -p '\tPress Enter to continue...' var");
#elif defined _WIN32
	system("pause");
#endif
}

//MAIN MENU
void mainMenu()
{
	clearScreen();

	int option = 0;

	//MAIN MENU
	cout << "\n";
	cout << "\t==============================================================\n";
	cout << "\t==                                                          ==\n";
	cout << "\t==                  LINEAR ALGEBRA SYSTEM                   ==\n";
	cout << "\t==                                                          ==\n";
	cout << "\t==              MS and PhD in Computer Science              ==\n";
	cout << "\t==               Federal University of Ceara                ==\n";
	cout << "\t==                       (MDCC/UFC)                         ==\n";
	cout << "\t==                                                          ==\n";
	cout << "\t==============================================================\n";
	cout << "\n";
	cout << "\tMAIN MENU\n";
	cout << "\n";
	cout << "\t1. Multiplying a Matrix by a Vector\n";
	cout << "\t2. Multiplying a Matrix by a Matrix\n";
	cout << "\t3. Gaussian Elimination (Partial Pivoting)\n";
	cout << "\t4. Gaussian Elimination (Full Pivoting)\n";
	cout << "\t5. Gauss-Jordan Elimination\n";
	cout << "\t6. LU Decomposition\n";
	cout << "\t7. Cholesky Decomposition\n";
	cout << "\t8. Reduced Row Echelon Form\n";
	cout << "\t9. The Gram-Schmidt Process\n";
	cout << "\t10. QR Decomposition (by Gram-Schmidt process)\n";
	cout << "\t11. QR Decomposition (by Reflection matrix)\n";
	cout << "\t12. QR Decomposition (by Rotation matrix)\n";
	cout << "\t13. Eigen-value and Eigen-vector (by Power Method)\n";
	cout << "\t14. Eigen-value and Eigen-vector (by Inverse Power Method)\n";
	cout << "\t15. Eigen-value and Eigen-vector (by Inverse Power Method with Displacement)\n";
	cout << "\t16. Eigen-value and Eigen-vector (by Householder Method)\n";
	cout << "\t0. EXIT\n";
	cout << "\n";
	cout << "\tEnter your choice: "; cin >> option;
	cout << endl;

	switch (option)
	{
	case 1:
		multiMatrixVector();
		break;
	case 2:
		multiMatrixMatrix();
		break;
	case 3:
		menuGaussPP();
		break;
	case 4:
		menuGaussFP();
		break;
	case 5:
		menuGaussJordan();
		break;
	case 6:
		menuLU();
		break;
	case 7:
		menuCholesky();
		break;
	case 8:
		menuRREF();
		break;
	case 9:
		menuOrth();
		break;
	case 10:
		menuQRbyGS();
		break;
	case 11:
		menuQRbyRef();
		break;
	case 12:
		menuQRbyRot();
		break;
	case 13:
		menuPowerMethod();
		break;
	case 14:
		menuInvPowerMethod();
		break;
	case 15:
		menuInvPowMethodDispl();
		break;
	case 16:
		menuHouseholderMethod();
		break;
	case 0:
		exit(0);
		break;
	default:
		cerr << "\n\t! INVALID OPTION !" << endl;
		pauseScreen();
		break;
	}
}

//Read Input File
bool readFile(Matrix& A, Vector& b)
{
	//Get the name of the input file
	string fName;
	//char fName[20];
	cout << "\n\tEnter the input file name [.dat]: "; cin >> fName;
	fName += ".dat";

	ifstream inFile(fName);
	if (!inFile.is_open())
	{
		cerr << "\n\t!!! FILE COULD NOT BE OPENED !!!\n" << endl;
		return false;
	}

	//Read a file line by line
	string line;
	bool sM = false;
	bool M = false;
	bool sV = false;
	bool V = false;

	int i = 0;
	while (!inFile.eof())
	{
		getline(inFile, line);

		/*
		 * ------------ %SIZE.M -------------
		 */

		if (line.substr(0, 7) == "%SIZE.M")
		{
			sM = true;
			//cout << "\n\tReading Size.M ..." << endl;
		}
		else {
			if (sM == true)
			{
				if (line.substr(0, 1) != "%")
				{
					istringstream s(line);
					s >> A.sizeR;
					s >> A.sizeC;
				}
				else sM = false;
			}
		}

		/*
		 * ------------ %MATRIX -------------
		 */

		if (line.substr(0, 7) == "%MATRIX")
		{
			i = 0;
			M = true;
			//cout << "\n\tReading Matrix ..." << endl;
		}
		else {
			if (M == true)
			{
				if (line.substr(0, 1) != "%")
				{
					istringstream s(line);
					for (int j = 0; j < A.sizeC; j++)
						s >> A.e[i][j];
					i++;
				}
				else M = false;
			}
		}

		/*
		 * ------------ %SIZE.V -------------
		 */

		if (line.substr(0, 7) == "%SIZE.V")
		{
			sV = true;
			//cout << "\n\tReading Size.V ..." << endl;
		}
		else {
			if (sV == true)
			{
				if (line.substr(0, 1) != "%")
				{
					istringstream s(line);
					s >> b.size;
				}
				else sV = false;
			}
		}

		/*
		 * ------------ %VECTOR -------------
		 */

		if (line.substr(0, 7) == "%VECTOR")
		{
			i = 0;
			V = true;
			//cout << "\n\tReading Vector ..." << endl;
		}
		else {
			if (V == true)
			{
				if (line.substr(0, 1) != "%")
				{
					istringstream s(line);
					s >> b.e[i];
					i++;
				}
				else V = false;
			}
		}

		//END OF FILE
	}
	return true;
}

//1. Multiplying a Matrix by a Vector
void multiMatrixVector()
{
	clearScreen();

	cout << "\n";
	cout << "\t==============================================================\n";
	cout << "\t==             MULTIPLYING A MATRIX BY A VECTOR             ==\n";
	cout << "\t==============================================================\n";
	cout << "\n";

	Matrix A;
	Vector x;
	if (readFile(A, x) == true)
	{
		cout << "\n\t*** MATRIX ***\n" << endl;
		A.print();
		cout << "\n\t*** VECTOR ***\n" << endl;
		x.print();

		Vector b = A.vectorMulti(x);
		cout << "\n\t*** RESULT ***\n" << endl;
		b.print();
	}

	pauseScreen();
}

//2. Multiplying a Matrix by a Matrix
void multiMatrixMatrix()
{
	clearScreen();

	cout << "\n";
	cout << "\t==============================================================\n";
	cout << "\t==             MULTIPLYING A MATRIX BY A MATRIX             ==\n";
	cout << "\t==============================================================\n";
	cout << "\n";

	//Get the name of the input file
	//string fName;
	char fName[20];

	cout << "\n\tEnter the input file name: "; cin >> fName;

	ifstream inFile(fName);
	if (!inFile.is_open())
	{
		cerr << "\n\t!!! FILE COULD NOT BE OPENED !!!\n" << endl;
		exit(EXIT_FAILURE);
	}

	Matrix L; //Matrix Left
	Matrix R; //Matrix Right
	//Read a file line by line
	string line;

	bool sL = false;
	bool mL = false;
	bool sR = false;
	bool mR = false;

	int i = 0;
	while (!inFile.eof())
	{
		getline(inFile, line);

		/*
		 * ------------ %SIZE.LEFT -------------
		 */

		if (line.substr(0, 7) == "%SIZE.LEFT")
		{
			sL = true;
			//cout << "\n\tReading Size.LEFT ..." << endl;
		}
		else {
			if (sL == true)
			{
				if (line.substr(0, 1) != "%")
				{
					istringstream s(line);
					s >> L.sizeR;
					s >> L.sizeC;
				}
				else sL = false;
			}
		}

		/*
		 * ------------ %MATRIX.LEFT -------------
		 */

		if (line.substr(0, 7) == "%MATRIX.LEFT")
		{
			i = 0;
			mL = true;
			//cout << "\n\tReading Matrix LEFT ..." << endl;
		}
		else {
			if (mL == true)
			{
				if (line.substr(0, 1) != "%")
				{
					istringstream s(line);
					for (int j = 0; j < L.sizeC; j++)
						s >> L.e[i][j];
					i++;
				}
				else mL = false;
			}
		}

		/*
		 * ------------ %SIZE.RIGHT -------------
		 */

		if (line.substr(0, 7) == "%SIZE.RIGHT")
		{
			sR = true;
			//cout << "\n\tReading Size.RIGHT ..." << endl;
		}
		else {
			if (sR == true)
			{
				if (line.substr(0, 1) != "%")
				{
					istringstream s(line);
					s >> R.sizeR;
					s >> R.sizeC;
				}
				else sR = false;
			}
		}

		/*
		 * ------------ %MATRIX.RIGHT -------------
		 */

		if (line.substr(0, 7) == "%MATRIX.RIGHT")
		{
			i = 0;
			mR = true;
			//cout << "\n\tReading Matrix RIGHT ..." << endl;
		}
		else {
			if (mR == true)
			{
				if (line.substr(0, 1) != "%")
				{
					istringstream s(line);
					for (int j = 0; j < R.sizeC; j++)
						s >> R.e[i][j];
					i++;
				}
				else mR = false;
			}
		}

		//END OF FILE
	}

	cout << "\n\t*** MATRIX LEFT ***\n" << endl;
	L.print();
	cout << "\n\t*** MATRIX RIGHT ***\n" << endl;
	R.print();

	Matrix Result = L.matrixMulti(R);
	cout << "\n\t*** RESULT ***\n" << endl;
	Result.print();

	pauseScreen();
}

//3. Gaussian Elimination (Partial Pivoting)
void menuGaussPP()
{
	clearScreen();

	cout << "\n";
	cout << "\t==============================================================\n";
	cout << "\t==                   GAUSSIAN ELIMINATION                   ==\n";
	cout << "\t==                     Partial Pivoting                     ==\n";
	cout << "\t==============================================================\n";
	cout << "\n";

	Matrix A;
	Vector b;

	if (readFile(A, b) == true)
		gaussPP(A, b);

	pauseScreen();
}

//4. Gaussian Elimination (Full Pivoting)
void menuGaussFP()
{
	clearScreen();

	cout << "\n";
	cout << "\t==============================================================\n";
	cout << "\t==                   GAUSSIAN ELIMINATION                   ==\n";
	cout << "\t==                       Full Pivoting                      ==\n";
	cout << "\t==============================================================\n";
	cout << "\n";

	Matrix A;
	Vector b;

	if (readFile(A, b) == true)
		gaussFP(A, b);

	pauseScreen();
}

//5. Gauss-Jordan Elimination
void menuGaussJordan()
{
	clearScreen();

	cout << "\n";
	cout << "\t==============================================================\n";
	cout << "\t==                 GAUSS-JORDAN ELIMINATION                 ==\n";
	cout << "\t==============================================================\n";
	cout << "\n";

	Matrix A;
	Vector b;

	if (readFile(A, b) == true)
		gaussJordan(A, b);

	pauseScreen();
}

//6. LU Decomposition
void menuLU()
{
	clearScreen();

	cout << "\n";
	cout << "\t==============================================================\n";
	cout << "\t==                     LU Decomposition                     ==\n";
	cout << "\t==============================================================\n";
	cout << "\n";

	Matrix A;
	Vector b;

	if (readFile(A, b) == true)
		LU(A, b);

	pauseScreen();
}

//7. Cholesky Decomposition
void menuCholesky()
{
	clearScreen();

	cout << "\n";
	cout << "\t==============================================================\n";
	cout << "\t==                  Cholesky Decomposition                  ==\n";
	cout << "\t==============================================================\n";
	cout << "\n";

	Matrix A;
	Vector b;

	if (readFile(A, b) == true)
		Cholesky(A, b);

	pauseScreen();
}

//8. Reduced Row Echelon Form
void menuRREF()
{
	clearScreen();

	cout << "\n";
	cout << "\t==============================================================\n";
	cout << "\t==                 Reduced Row Echelon Form                 ==\n";
	cout << "\t==============================================================\n";
	cout << "\n";

	Matrix A;
	Vector b;

	if (readFile(A, b) == true)
	{
		Matrix R;
		R = RREF(A);
	}

	pauseScreen();
}

//9. The Gram-Schmidt Process
void menuOrth()
{
	clearScreen();

	cout << "\n";
	cout << "\t==============================================================\n";
	cout << "\t==                ORTHONORMAL BASIS VECTORS                 ==\n";
	cout << "\t==                (by GRAM-SCHMIDT Process)                 ==\n";
	cout << "\t==============================================================\n";
	cout << "\n";

	Matrix A;
	Vector b;

	if (readFile(A, b) == true)
		orthonormal(A);

	pauseScreen();
}

//10. QR Decomposition (by Gram-Schmidt process)
void menuQRbyGS()
{
	clearScreen();

	cout << "\n";
	cout << "\t==============================================================\n";
	cout << "\t==                    QR Decomposition                      ==\n";
	cout << "\t==                (by GRAM-SCHMIDT Process)                 ==\n";
	cout << "\t==============================================================\n";
	cout << "\n";

	Matrix A;
	Vector b;

	if (readFile(A, b) == true)
		QRbyGS(A, b);

	pauseScreen();
}

//11. QR Decomposition (by Reflection matrix)
void menuQRbyRef()
{
	clearScreen();

	cout << "\n";
	cout << "\t==============================================================\n";
	cout << "\t==                    QR Decomposition                      ==\n";
	cout << "\t==                 (by Reflection Matrix)                   ==\n";
	cout << "\t==============================================================\n";
	cout << "\n";

	Matrix A;
	Vector b;

	if (readFile(A, b) == true)
		QRbyRef(A, b);

	pauseScreen();
}

//12. QR Decomposition (by Rotation matrix)
void menuQRbyRot()
{
	clearScreen();

	cout << "\n";
	cout << "\t==============================================================\n";
	cout << "\t==                    QR Decomposition                      ==\n";
	cout << "\t==                  (by Rotation Matrix)                    ==\n";
	cout << "\t==============================================================\n";
	cout << "\n";

	Matrix A;
	Vector b;

	if (readFile(A, b) == true)
		QRbyRot(A, b);

	pauseScreen();
}

//13. Eigenvalue (by Power Method)
void menuPowerMethod()
{
	clearScreen();

	cout << "\n";
	cout << "\t==============================================================\n";
	cout << "\t==               Eigenvalue and Eigenvector                 ==\n";
	cout << "\t==                    (by Power Method)                     ==\n";
	cout << "\t==============================================================\n";
	cout << "\n";

	Matrix A;
	Vector b;

	if (readFile(A, b) == true)
		powerMethod(A, b);

	pauseScreen();
}

//14. Eigenvalue (by Inverse Power Method)
void menuInvPowerMethod()
{
	clearScreen();

	cout << "\n";
	cout << "\t==============================================================\n";
	cout << "\t==               Eigenvalue and Eigenvector                 ==\n";
	cout << "\t==                (by Inverse Power Method)                 ==\n";
	cout << "\t==============================================================\n";
	cout << "\n";

	Matrix A;
	Vector b;

	if (readFile(A, b) == true)
		invPowerMethod(A, b);

	pauseScreen();
}

//15. Eigen-value and Eigen-vector (by Inverse Power Method with Displacement)
void menuInvPowMethodDispl()
{
	clearScreen();

	cout << "\n";
	cout << "\t==============================================================\n";
	cout << "\t==               Eigenvalue and Eigenvector                 ==\n";
	cout << "\t==       (by Inverse Power Method with Displacement)        ==\n";
	cout << "\t==============================================================\n";
	cout << "\n";

	Matrix A;
	Vector b;

	if (readFile(A, b) == true)
		invPowMethodDispl(A, b);

	pauseScreen();
}

//16. Eigen-value and Eigen-vector (by Householder Method)
void menuHouseholderMethod()
{
	clearScreen();

	cout << "\n";
	cout << "\t==============================================================\n";
	cout << "\t==               Eigenvalue and Eigenvector                 ==\n";
	cout << "\t==                 (by Householder Method)                  ==\n";
	cout << "\t==============================================================\n";
	cout << "\n";

	Matrix A;
	Vector b;

	if (readFile(A, b) == true)
		HouseholderMethod(A);

	pauseScreen();
}
