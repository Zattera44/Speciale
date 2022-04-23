#pragma once
#include <vector>

double trapezoid(double (*func)(double), double a, double b, int n);
double hyperGeom(double a, double b, double c, double z);
double G(double gamma, double t, double u);


std::vector<std::vector<double>> cholesky(const std::vector<std::vector<double>>& A);
std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>>& A);
std::vector<std::vector<double>> matrixMult(const std::vector<std::vector<double>>& A
										   ,const std::vector<std::vector<double>>& B);
std::vector<double> vectorMult(const std::vector<std::vector<double>>& A, const std::vector<double>& a);



void printMatrix(const std::vector<std::vector<double>>& M);
void printVector(const std::vector<double>& m);


double normCDF(double x);
double normPDF(double x);


