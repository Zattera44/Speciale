#include "Simulation.h"
#include <random>

std::vector<std::vector<double>> covMatrix(const double& H, const double& T, const int& n)
{
	double dt;
	std::vector<std::vector<double>> cov(n, std::vector<double>(n,0));

	dt = (double) T / n;
	for (size_t i = 0; i < cov.size(); i++)
	{
		for (size_t j = 0; j < cov.size(); j++)
		{
			cov[i][j] = 0.5 * (pow((i+1) * dt, 2 * H) + pow((j+1) * dt, 2 * H) - pow(fabs((i+1) * dt - (j+1) * dt), 2 * H));
		}

	}
	return cov;
}



std::vector<double> BSTest(double spot_,double strike_, double r_,double sigma_, double T_, double rand)
{
	std::vector<double> result(5);

	Number spot(spot_);
	Number strike(strike_);
	Number r(r_);
	Number sigma(sigma_);
	Number T(T_);


	Number price(exp(-1*r*T) * max(spot * exp( (r - 0.5 * sigma*sigma ) * T + sigma * sqrt(T) * rand) - strike,0.0));


	price.propogateAdjoints();


	result[0] = price.getNode()->result();
	result[1] = spot.adjoint();
	result[2] = r.adjoint();
	result[3] = sigma.adjoint();
	result[4] = T.adjoint();	

	sigma.resetTape();

	return result;
}


std::vector<double> BSSim(double spot, double strike, double r, double sigma, double T, int N)
{
	std::vector<double> temp(5);
	std::vector<double> results(5);

	results[0] = 0;
	results[1] = 0;
	results[2] = 0;
	results[3] = 0;
	results[4] = 0;

	double random;

	std::mt19937_64 rand(30);

	std::normal_distribution<double> norm(0.0, 1.0);


	for (int i = 0; i < N; i++)
	{
		random = norm(rand);
		temp = BSTest(spot, strike, r, sigma, T, random);
	
		results[0] += temp[0];
		results[1] += temp[1];
		results[2] += temp[2];
		results[3] += temp[3];
		results[4] += temp[4];


	}

	results[0] = (double) results[0] / N;
	results[1] = (double) results[1] / N;
	results[2] = (double) results[2] / N;
	results[3] = (double) results[3] / N;
	results[4] = (double) results[4] / N;

	return results;
}



std::vector<double> BSFormula(double spot_, double strike_, double r_, double sigma_, double T_)
{
	std::vector<double> result(5);
	
	Number spot(spot_);
	Number strike(strike_);
	Number r(r_);
	Number sigma(sigma_);
	Number T(T_);


	Number d1((log(spot / strike) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T)));
	Number d2(d1 - sigma * sqrt(T));

	Number price(spot * CDF(d1) - exp(-1 * r * T) * strike * CDF(d2));

	price.propogateAdjoints();

	result[0] = price.getNode()->result();
	result[1] = spot.adjoint();
	result[2] = sigma.adjoint();
	result[3] = r.adjoint();
	result[4] = T.adjoint();

	return result;
}
