#include <iostream>
#include <random>
#include <vector>

#include "Utils.h"
#include "Simulation.h"




int main()
{


	double spot = 100;
	double strike = 90;
	double r = 0;
	double sigma = 0.15;
	double T = 1;
	int N = 10000000;


	std::vector<double> results(5);


	//results = BSSim(spot, strike, r, sigma, T, N);

	results = BSFormula(spot, strike, r, sigma, T);

	for (int i = 0; i < results.size(); i++)
	{
		std::cout << results[i] << std::endl;
	}


}
