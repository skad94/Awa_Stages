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

unique_ptr<vector<double>> StandardBM_QuickMethod(int numTimeSteps, double maturity);

double BlackScholesCallPrice(double maturity, double strike, double volatility, 
	double shortRate, double S0);


// rHeston

double DiffusionRoughHeston(int numTimeSteps, double maturity, double P0, 
	double Y0, double alpha, double lambdaStar, double mu, double beta, double C, 
	double phi1Norm, double phi2Norm);

double DiffusionRoughHeston_AbiJaber(int numTimeSteps, double maturity, double S0,
	double V0, double H, double lambda, double theta, double nu, double rho);

unique_ptr<vector<double>> VariancePath_RoughHeston_AbiJaber(int numTimeSteps, double maturity, 
	double V0, double H, double lambda, double theta, double nu);

unique_ptr<vector<double>> CalibrateRoughHeston_AbiJaber(double VIX, 
	const unique_ptr<vector<double>>& marketPrices, const unique_ptr<vector<double>>& maturities, 
	const unique_ptr<vector<double>>& strikes, string financialProduct, int numMonteCarloSimulations, 
	int numTimeSteps, double shortRate);


// rBergomi

double IntegrandBergomi(double x, double s, double gamma);

double IntegrateCovBergomi(double(*f)(double, double, double), double x,
	double gamma, double lowerBound, double upperBound);

double Covariance_Wtilde(double s, double t, double HurstExponent);

double Covariance_StandardBM(double s, double t);

double Covariance_Wtilde_Z(double s, double t, double HurstExponent, double rho);

unique_ptr<vector<double>> Trajectory_Wtilde_and_Z(int numTimeSteps, const unique_ptr<cSquareMatrix>& sigma);

double DiffusionRoughBergomi(int n, double maturity, double S0, double v0,
	double eta, const unique_ptr<cSquareMatrix>& sigma);

double MonteCarlo_rBergomi_Call(int numMonteCarloSimulations, int numTimeSteps, double maturity,
	double S0, double v0, double eta, const unique_ptr<cSquareMatrix>& sigma, 
	double strike, double shortRate);


// Lifted Heston

double Weight_LiftedHeston(int n, int i, double rn, double H);

double ReversionSpeed_LiftedHeston(int n, int i, double rn, double H);

double InitialForwardVarianceCurve(double t, double V0, double theta, int n, double rn, double H);

double DiffusionLiftedHeston(int numTimeSteps, double maturity, int n, double rn, double S0,
	double V0, double theta, double rho, double nu, double lambda, double H);

unique_ptr<vector<double>> CalibrateLiftedHeston(const unique_ptr<vector<double>>& marketPrices,
	const unique_ptr<vector<double>>& maturities, const unique_ptr<vector<double>>& strikes,
	string financialProduct, int numMonteCarloSimulations, int numTimeSteps, double shortRate);

double ImplicitVolLiftedHeston(double LiftedHestonPrice, double maturity, double strike, 
	double shortRate);

unique_ptr<vector<double>> GenerateVolSurf_LiftedHeston(double shortRate);


// Monte Carlo pricing

double MonteCarlo_pricing(string financialProduct, string model, int numMonteCarloSimulations, 
	int numTimeSteps, double maturity, double strike, double shortRate, int n = 20, 
	double rn = 2.5, const unique_ptr<cSquareMatrix>& sigma = unique_ptr<cSquareMatrix>(new cSquareMatrix(0)));

double VarianceDerivatives_pricing_rHeston(string financialProduct, int numMonteCarloSimulations, 
	int numTimeSteps, double maturity, double V0, double H, double lambda, double theta, double nu,
	double strike, double shortRate);

double VIX_options_pricing_rHeston(string financialProduct, int numMonteCarloSimulations,
	int numTimeSteps, double maturity, double V0, double H, double lambda, double theta, double nu,
	double strike, double shortRate);

#endif