#pragma once
#include <vector>

double trapezoid(double (*func)(double), double a, double b, int n);
std::vector<std::vector<double>> cholesky(std::vector<std::vector<double>> A);
std::vector<std::vector<double>> transpose(std::vector<std::vector<double>> A);
std::vector<std::vector<double>> squareMult(std::vector<std::vector<double>> A
										   ,std::vector<std::vector<double>> B);


void printMatrix(std::vector<std::vector<double>> M);