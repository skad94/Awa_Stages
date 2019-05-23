#ifndef FUNCTIONS_HH
#define FUNCTIONS_HH

#include "cSquareMatrix.hh"

double Integrate(double(*f)(double, double, int, double), double H,
	double T, int k, double a, double b);

double CovarianceFBM(double s, double t, double H);

unique_ptr<vector<double>> SamplePathFBM_Cholesky(double maturity, int n, double H);

double CosinusKL(double t, int k, double T);

double SinusKL(double t, int k, double T);

double F_KL(double H, double T, int k, double u);

double G_KL(double H, double T, int k, double u);

double Coeff_C0_KL(double H, double T);

double Coeff_Ck_KL(double H, double T, int k);

unique_ptr<vector<double>> SamplePathFBM_KL(double maturity, int n, double H, int N);

#endif