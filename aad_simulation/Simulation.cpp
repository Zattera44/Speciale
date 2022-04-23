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


double BSCall(double spot, double strike, double r, double sigma, double T)
{
	std::mt19937_64 rand(30);

	std::normal_distribution<double> norm(0.0, 1.0);

	double rng;
	double sum = 0;
	int n = 100000;


	for (int i = 0; i < n; i++)
	{
		rng = norm(rand);
		sum = sum + fmax(spot * exp((r - 0.5 * sigma * sigma) * T + sigma * pow(T,0.5) * rng) - strike, 0.0);
	}

	return exp(-r * T) * sum / n;
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



std::vector<double> rBergomiCall(double H, double xi, double eta, double spot, double v0, double T, int numSteps)
{
	std::vector<double> res(5);
	//Generer Z og Wtilde

	//Alt herunder skal med i AAD
		//Simuler V

		//Simuler S



		//Beregn payoff og greeks
	return res;
}



std::vector < std::vector<double> > rBergomiCov(double H,double rho, double T, int numSteps)
{
	
	std::vector<std::vector<double>> covMatrix(2*numSteps, std::vector<double>(2*numSteps,0));
	double dt = (double)T / numSteps;


	//Dan kovariansmatrix, evt. blok for blok

	//Kovariansmatrix for Z
	for (int i = 0; i < numSteps; i++)
	{
		for (int j = 0; j < numSteps; j++) 
		{
			covMatrix[i][j] = fmax(i * dt, j * dt);
		}
	}


	//Kovariansmatrix for Wtilde
	for (int i = numSteps; i < 2*numSteps; i++)
	{
		for (int j = numSteps; j < 2*numSteps; j++)
		{
			if (i < j)
			{
				covMatrix[i][j] = G(H - 0.5, 0.5*i * dt, 0.5*j * dt);
			}
			if (j > i)
			{
				covMatrix[i][j] = G(0.5 - H, 0.5*j * dt, 0.5*i * dt);
			}
			else
			{
				covMatrix[i][j] = pow(0.5*i * dt, 2 * H);
			}
		}
	}

	//øvre kovariansmatrix mellem Wtilde og Z

	double Dh = sqrt(2 * H) / (H + 0.5);

	for (int i = 0; i < numSteps; i++)
	{
		for (int j = numSteps; j < 2*numSteps; j++)
		{

		}
	}




	return covMatrix;
}