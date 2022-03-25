#include "Simulation.h"
#include <random>
#include "Node.h"


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



//std::vector<double> BStest(Number spot, Number strike, Number r, Number sigma, Number T, double rand)
//{
//	std::vector<double> result(5);
//
//
//	Number St(spot * exp((r - 0.5 * square(sigma)) * T + sigma * sqrt(T) * rand));
//
//	spot.propogateAdjoints();
//
//
//	result[0] = St.getNode()->result();
//	result[1] = spot.adjoint();
//	result[2] = r.adjoint();
//	result[3] = sigma.adjoint();
//	result[4] = T.adjoint();
//
//
//	return result;
//}


//std::vector<double> BSSim(Number spot, Number strike, Number r, Number sigma, Number T, int N)
//{
//	std::vector<double> results(5);
//
//	for (int i = 0; i < N; i++)
//	{
//
//		results = BSTest(spot, strike, r, sigma, T, 10);
//
//		
//
//
//	}
//
//}