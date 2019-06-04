#ifndef FUNCTIONS_HH
#define FUNCTIONS_HH

#include "cSquareMatrix.hh"
#include <cmath>
#include <random>
#include <ctime>
#include <fstream>

double IntegrateForSeriesExpansionMethod(double(*f)(double, double, int, double), double HurstExponent,
	double maturityKL, int k_index_KL, double lowerBound, double upperBound);

double CovarianceFBM(double s, double t, double H);

unique_ptr<vector<double>> SamplePathFBM_Cholesky(double maturity, int numTimeSteps, double H);

double CosinusKL(double t, int k_index_KL, double maturityKL);

double SinusKL(double t, int k_index_KL, double maturityKL);

double F_KL(double H, double maturityKL, int k_index_KL, double variableOfIntegration);

double G_KL(double H, double maturityKL, int k_index_KL, double variableOfIntegration);

double Coeff_C0_KL(double H, double maturityKL);

double Coeff_Ck_KL(double H, double maturityKL, int k_index_KL);

unique_ptr<vector<double>> SamplePathFBM_KL(double maturity, int numTimeSteps, double H, int upperBoundInSumKL);

void ExportData(const unique_ptr<vector<double>>& data, string fileName);

unique_ptr<vector<double>> StandardBM(int numTimeSteps, const unique_ptr<cSquareMatrix>& sigma);


// rHeston

double DiffusionRoughHeston(int n, double maturity, double P0, double Y0, double alpha, 
	double lambdaStar, double mu, double beta, double C, double phi1Norm, double phi2Norm, 
	const unique_ptr<cSquareMatrix>& sigma);


// rBergomi

double IntegrandBergomi(double x, double s, double gamma);

double IntegrateCovBergomi(double(*f)(double, double, double), double x,
	double gamma, double a, double b);

double Covariance_Wtilde(double u, double v, double H);

double Covariance_StandardBM(double u, double v);

double Covariance_Wtilde_Z(double u, double v, double H, double rho);

unique_ptr<vector<double>> Trajectory_Wtilde_and_Z(int n, const unique_ptr<cSquareMatrix>& sigma);

double DiffusionRoughBergomi(int n, double maturity, double S0, double v0, double eta,
	const unique_ptr<cSquareMatrix>& sigma);


#endif