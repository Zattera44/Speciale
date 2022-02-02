#include "Utils.h"



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