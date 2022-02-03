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
    std::vector<std::vector<double>> A = { {2 , 1 ,1}, {1,2,1},{1,1,2}};
    std::vector<std::vector<double>> B(3, std::vector<double>(3, 1));
    
    std::cout << "A" << std::endl;
    
    printMatrix(A);
    

    std::cout << "L" << std::endl;

    B = cholesky(A);

    printMatrix(B);


    std::cout << std::endl << "L L^T" << std::endl;;

    printMatrix(squareMult(B,transpose(B)));
}
