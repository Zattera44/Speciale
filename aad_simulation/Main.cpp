#include <iostream>
#include <vector>
#include <chrono>

#include "Utils.h"
#include "Simulation.h"


int main()
{
	double spot = 100;
	double strike = 90;
	double r = 0.03;
	double T = 1;
	double sigma = 0.2;
	int N = 100000;

	std::vector<double> res;


	auto start = std::chrono::steady_clock::now();

	res = BSSim(spot, strike, r, sigma, T, N);

	auto end = std::chrono::steady_clock::now();

	std::cout << "Time spent: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()  << std::endl;

	for (auto i : res)
	{
		std::cout << i << std::endl;
	}


	start = std::chrono::steady_clock::now();

	double h = 0.001;

	std::cout << "Finite difference" << std::endl;

	std::cout << BSCall(spot, strike, r, sigma, T) << std::endl;
	std::cout << (BSCall(spot + h, strike, r, sigma, T) - BSCall(spot - h, strike, r, sigma, T)) / (2*h) << std::endl;
	std::cout << (BSCall(spot, strike, r + h, sigma, T) - BSCall(spot, strike, r - h, sigma, T)) / (2*h) << std::endl;
	std::cout << (BSCall(spot, strike, r, sigma + h, T) - BSCall(spot, strike, r, sigma - h, T)) / (2*h) << std::endl;
	std::cout << (BSCall(spot, strike, r, sigma, T + h) - BSCall(spot, strike, r, sigma, T - h)) / (2*h) << std::endl;

	end = std::chrono::steady_clock::now();


	std::cout << "Time spent: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << std::endl;


}
