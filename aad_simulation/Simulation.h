#pragma once
#include <vector>
#include <random>
#include "Number.h"

std::vector<std::vector<double>> covMatrix(const double& H, const double& T, const int& n);


std::vector<double> BSTest(double spot, double strike,double r, double sigma, double T, double rand);

std::vector<double> BSSim(double spot, double strike, double r, double sigma, double T, int N);


std::vector<double> BSFormula(double spot, double strike, double r, double sigma, double T);
