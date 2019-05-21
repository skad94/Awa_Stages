#include "Functions.hh"
#include <cmath>

double 
CovarianceFMB(double s, double t, double H)
{
	return pow(s, 2 * H) + pow(t, 2 * H) - pow(abs(t - s), 2 * H);
}

