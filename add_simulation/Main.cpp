// SpecialeSimulation.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "Utils.h"

double f(double x)
{
    return pow(x, 2);
}


int main()
{
    double a, b, n;

    a = 0;
    b = 4;
    n = 1000;

    std::cout << trapezoid(f, a, b, n) << std::endl;

    std::cout << "Hello World!\n";
}
