#include "Utils.h"
#include <iostream>

void printMatrix(std::vector<std::vector<double>> M)
{
	for (int i = 0; i < M.size(); i++)
	{
		for (int j = 0; j < M[i].size(); j++)
		{
			std::cout << M[i][j];
		
			if (j == M[i].size() - 1 ) { std::cout << std::endl; }
			else { std::cout << " "; }
		}
		
	}
}


// Numerical integration using the trapezoid rule
double trapezoid(double (*func)(double), double a, double b, int n)
{
	double x;
	double sum;
	double dx;


	dx = (double) (b - a) / n;
	x = a + dx;

	sum = (func(b) + func(a)) / 2;

	for (int i = 1; i < n; i++)
	{
		sum += func(x); 
		x += dx;
	}

	sum = sum * dx;

	return sum;
}



//Returns the transpose of a matrix A
std::vector<std::vector<double>> transpose(std::vector<std::vector<double>> A)
{
	std::vector<std::vector<double>> ATrans(A.size(), std::vector<double>(A[0].size(), 0));

	for (int i = 0; i < A.size(); i++)
	{
		for (int j = 0; j < A[0].size(); j++)
		{
			ATrans[i][j] = A[j][i];
		}
	}
	return ATrans;
}



//Returns the product of two matrices A and B
std::vector<std::vector<double>> squareMult(std::vector<std::vector<double>> A
	, std::vector<std::vector<double>> B)
{
	std::vector<std::vector<double>> C(A.size(), std::vector<double>(A[0].size(), 0));


	for (int i = 0; i < C.size(); i++)
	{
		for (int j = 0; j < C[i].size(); j++)
		{
			double sum = 0;
			for (int k = 0; k < A[i].size(); k++)
			{
				sum += A[i][k] * B[k][j];
			}
			C[i][j] = sum;
		}
	}
	return C;
}



//Returns the Cholesky decomposition of a positive definite matrix
std::vector<std::vector<double>> cholesky(std::vector<std::vector<double>> A)
{
	std::vector<std::vector<double>> L(A.size(), std::vector<double>(A[0].size(),0));
	for (int i = 0; i < L.size(); i++)
	{
		for (int j = 0; j <= i ; j++)
		{
			double sum = 0;
			if (i == j)
			{
				for (int k = 0; k < j; k++)
				{
					sum += pow(L[j][k], 2);
				}
				L[i][j] = pow(A[i][j] - sum,0.5);
			}
			else
			{
				for (int k = 0; k < j; k++)
				{
					sum += L[i][k] * L[j][k];
				}
				L[i][j] = 1 / L[j][j] * (A[i][j] - sum);
			}
		}
	}
	return L;
}