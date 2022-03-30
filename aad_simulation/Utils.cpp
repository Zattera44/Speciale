#include "Utils.h"
#include <iostream>

void printMatrix(const std::vector<std::vector<double>>& M)
{
	for (size_t i = 0; i < M.size(); i++)
	{
		for (size_t j = 0; j < M[i].size(); j++)
		{
			std::cout << M[i][j];
		
			if (j == M[i].size() - 1 ) { std::cout << std::endl; }
			else { std::cout << " "; }
		}
		
	}
}


void printVector(const std::vector<double>& m)
{
	std::cout << std::endl;
	for (size_t i = 0; i < m.size(); i++)
	{
		std::cout << m[i] << std::endl;
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
std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>>& A)
{
	std::vector<std::vector<double>> ATrans(A.size(), std::vector<double>(A[0].size(), 0));

	for (size_t i = 0; i < A.size(); i++)
	{
		for (size_t j = 0; j < A[0].size(); j++)
		{
			ATrans[i][j] = A[j][i];
		}
	}
	return ATrans;
}



//Returns the product of two matrices A and B
std::vector<std::vector<double>> matrixMult(const std::vector<std::vector<double>>& A
	, const std::vector<std::vector<double>>& B)
{
	std::vector<std::vector<double>> C(A.size(), std::vector<double>(B[0].size(), 0));


	for (size_t i = 0; i < C.size(); i++)
	{
		for (size_t j = 0; j < C[i].size(); j++)
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


//returns the result of A*a
std::vector<double> vectorMult(const std::vector<std::vector<double>>& A, const std::vector<double>& a)
{
	std::vector<double> res(A[0].size(), 0);
	double sum;


	for (size_t i = 0; i < res.size(); i++)
	{
		sum = 0;
		for (size_t k = 0; k < res.size(); k++)
		{
			sum += a[k] * A[i][k];
		}
		res[i] = sum;
		sum = 0;
	}
	return res;
}




//Returns the Cholesky decomposition of a positive definite matrix
std::vector<std::vector<double>> cholesky(const std::vector<std::vector<double>>& A)
{
	std::vector<std::vector<double>> L(A.size(), std::vector<double>(A[0].size(),0));
	for (size_t i = 0; i < L.size(); i++)
	{
		for (size_t j = 0; j <= i ; j++)
		{
			double sum = 0;
			if (i == j)
			{
				for (size_t k = 0; k < j; k++)
				{
					sum += pow(L[j][k], 2);
				}
				L[i][j] = pow(A[i][j] - sum,0.5);
			}
			else
			{
				for (size_t k = 0; k < j; k++)
				{
					sum += L[i][k] * L[j][k];
				}
				L[i][j] = 1 / L[j][j] * (A[i][j] - sum);
			}
		}
	}
	return L;
}




double normCDF(double x)
{
	static const double RT2PI = sqrt(4.0 * acos(0.0));

	static const double SPLIT = 7.07106781186547;

	static const double N0 = 220.206867912376;
	static const double N1 = 221.213596169931;
	static const double N2 = 112.079291497871;
	static const double N3 = 33.912866078383;
	static const double N4 = 6.37396220353165;
	static const double N5 = 0.700383064443688;
	static const double N6 = 3.52624965998911e-02;
	static const double M0 = 440.413735824752;
	static const double M1 = 793.826512519948;
	static const double M2 = 637.333633378831;
	static const double M3 = 296.564248779674;
	static const double M4 = 86.7807322029461;
	static const double M5 = 16.064177579207;
	static const double M6 = 1.75566716318264;
	static const double M7 = 8.83883476483184e-02;

	const double z = fabs(x);
	double c = 0.0;

	if (z <= 37.0)
	{
		const double e = exp(-z * z / 2.0);
		if (z < SPLIT)
		{
			const double n = (((((N6 * z + N5) * z + N4) * z + N3) * z + N2) * z + N1) * z + N0;
			const double d = ((((((M7 * z + M6) * z + M5) * z + M4) * z + M3) * z + M2) * z + M1) * z + M0;
			c = e * n / d;
		}
		else
		{
			const double f = z + 1.0 / (z + 2.0 / (z + 3.0 / (z + 4.0 / (z + 13.0 / 20.0))));
			c = e / (RT2PI * f);
		}
	}
	return x <= 0.0 ? c : 1 - c;
}


double normPDF(double x)
{
	return 1 / (sqrt(2 * 3.14159265358979323846264338327)) * exp(-0.5 * x * x);
}