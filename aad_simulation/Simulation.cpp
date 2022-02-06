#include "Simulation.h"


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
