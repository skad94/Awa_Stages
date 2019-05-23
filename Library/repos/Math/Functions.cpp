#define _USE_MATH_DEFINES 
#include "Functions.hh"
#include <cmath>
#include <random>
#include <ctime>

double 
Integrate(double(*f)(double, double, int, double), double H, double T, int k, double a, double b) 
{//Aproximation to an integral using rectangle method over [a, b], used to simulate a fractional 
 //brownian motion sample path by Series Expansion method
	int n = 10000;
	double dt = (b - a) / n;
	double res = 0;
	for (int i = 0; i < n; i++)
	{
		res += dt * f(H, T, k, a + (i + 1.) * dt);
	}
	return res;
}
	
double 
CovarianceFBM(double s, double t, double H)
{//Covariance function of a fractional brownian motion with Hurst exponent H
	return pow(s, 2 * H) + pow(t, 2 * H) - pow(abs(t - s), 2 * H);
}

unique_ptr<vector<double>> 
SamplePathFBM_Cholesky(double maturity, int n, double H)
{//Simulation of a fractional brownian motion sample path with exponent H by Cholesky decomposition
	double dt = maturity / n;
	default_random_engine generator;
	generator.seed(time(NULL));
	normal_distribution<double> distribution(0, 1);
	vector<double> gaussianSample(n, 0);
	cSquareMatrix covarianceMatrix(n);
	for (int i = 0; i < n; i++)
	{
		gaussianSample[i] = distribution(generator);
		for (int j = 0; j < n; j++)
			covarianceMatrix(i, j) = CovarianceFBM((i + 1.) * dt, (j + 1.) * dt, H);
	}
	unique_ptr<cSquareMatrix> sigma = covarianceMatrix.Cholesky();
	unique_ptr<vector<double>> res(*sigma * gaussianSample);
	return res;
}

double 
CosinusKL(double t, int k, double T)
{//This function will be useful to simulate a fractional brownian motion sample path by 
 //Series Expansion method
	return cos(M_PI * t * k / T);
}

double
SinusKL(double t, int k, double T)
{//This function will be useful to simulate a fractional brownian motion sample path by 
 //Series Expansion method
	return sin(M_PI * t * k / T);
}

double 
F_KL(double H, double T, int k, double u)
{//This function will be useful to simulate a fractional brownian motion sample path by 
 //Series Expansion method
	return pow(u, 2 * H)* CosinusKL(u, k, T);
}

double
G_KL(double H, double T, int k, double u)
{//This function will be useful to simulate a fractional brownian motion sample path by 
 //Series Expansion method
	return pow(u, 2 * H - 2) * CosinusKL(u, k, T);
}

double 
Coeff_C0_KL(double H, double T)
{//This function will be useful to simulate a fractional brownian motion sample path by 
 //Series Expansion method
	if (H < 0.5)
		return 0;
	return H * pow(T, 2 * H - 2);
}

double 
Coeff_Ck_KL(double H, double T, int k)
{//This function will be useful to simulate a fractional brownian motion sample path by 
 //Series Expansion method
	if (H < 0.5)
		return (2 / T) * Integrate(F_KL, H, T, k, 0, T);
	return (-4 * H * (2 * H - 1) * T / pow(k * M_PI, 2)) * Integrate(G_KL, H, T, k, 0, T);
}

unique_ptr<vector<double>> 
SamplePathFBM_KL(double maturity, int n, double H, int N)
{//Simulation of a fractional brownian motion sample path with Hurst exponent H by Series Expansion.
 //N is the number of terms in the infinite sum in Karhunen-Loève theorem.
	double dt = maturity / n;
	default_random_engine generator;
	generator.seed(time(NULL));
	normal_distribution<double> distribution(0, 1);
	vector<double> gaussianSample(2 * N + 1, 0);
	vector<double> coeff_Ck(N, 0);
	unique_ptr<vector<double>> FBM(new vector<double>(n, 0));
	for (int i = 0; i < N; i++) //We generate the gaussian sample and we compute the coefficients ck
	{
		coeff_Ck[i] = Coeff_Ck_KL(H, maturity, i + 1);
		gaussianSample[2 * i] = distribution(generator);
		gaussianSample[2 * i + 1] = distribution(generator);
	}
	for (int i = 0; i < n; i++)
	{
		double res = sqrt(Coeff_C0_KL(H, maturity)) * (i + 1.) * dt * gaussianSample[0];
		for (int k = 0; k < N; k++)
		{
			try //DEBUG
			{
				if (coeff_Ck[k] >= 0) throw "sqrt(negative number)";
			}
			catch (const char* err)
			{
				cerr << err << endl;
				exit(1);
			}
			double temp = sqrt(-coeff_Ck[k] / 2);
			res += temp * (SinusKL((i + 1.) * dt, k + 1, maturity) * gaussianSample[2 * k + 1]
				+ (1 - CosinusKL((i + 1.) * dt, k + 1, maturity)) * gaussianSample[2 * k + 2]);
		}
		(*FBM)[i] = res;
	}
	return FBM;
}