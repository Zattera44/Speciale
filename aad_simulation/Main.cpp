#include <iostream>
#include <random>
#include <vector>


#include "Utils.h"
#include "Simulation.h"


int main()
{

	double a, b, c, z;

	a = 1;
	b = 1;
	c = 2;
	z = 0.5;

	double gamma = 0.4;
	double x = 1.5;

	std::cout <<  G(gamma, 1,2) / (pow(2, 0.2) * pow(1, 0.2)) << std::endl;

}
