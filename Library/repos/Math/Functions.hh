#ifndef FUNCTIONS_HH
#define FUNCTIONS_HH

#include "cSquareMatrix.hh"
#include <cmath>
#include <random>
#include <ctime>
#include <fstream>

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

void ExportData(const unique_ptr<vector<double>>& data, string fileName);

unique_ptr<vector<double>> StandardBM(double maturity, int n, const unique_ptr<cSquareMatrix>& sigma);


// rHeston

double DiffusionRoughHeston(int n, double maturity, double P0, double Y0, double alpha, 
	double lambdaStar, double mu, double beta, double C, double phi1Norm, double phi2Norm, 
	const unique_ptr<cSquareMatrix>& sigma);

// rBergomi

double IntegrandBergomi(double x, double s, double gamma);

double IntegrateCovBergomi(double(*f)(double, double, double), double x,
	double gamma, double a, double b);

double Covariance_Wtilde(double u, double v, double H);

double Covariance_StandardMB(double u, double v);

double Covariance_Wtilde_Z(double u, double v, double H, double rho);



#endif