#include "Model.h"
#include <math.h>

double BlackScholes::increment(double spot, double dSpot)
{
	return spot * exp( (r - 0.5 * sigma)*dSpot + sigma * pow(dSpot,0.5) * 1 );
}