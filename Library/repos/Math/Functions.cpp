#define _USE_MATH_DEFINES 
#include "Functions.hh"
#include <iostream>
using namespace std;

double 
IntegrateForSeriesExpansionMethod(
	double(*f)(double, double, int, double), 
	double HurstExponent, 
	double maturityKL, 
	int k_index_KL, 
	double lowerBound, 
	double upperBound) 
{//Aproximation to an integral using rectangle method over [lowerBound, upperBound] used to 
 //simulate a fractional brownian motion sample path by Series Expansion method
	int numSteps = 50000;
	double dt = (upperBound - lowerBound) / numSteps;
	double res = 0;
	for (int i = 0; i < numSteps; i++)
	{
		res += dt * f(HurstExponent, maturityKL, k_index_KL, lowerBound + (i + 1/2) * dt);
	}
	return res;
}
	
double 
CovarianceFBM(double s, double t, double H)
{//Covariance function of a fractional brownian motion with Hurst exponent H
	return pow(s, 2 * H) + pow(t, 2 * H) - pow(abs(t - s), 2 * H);
}

unique_ptr<vector<double>> 
SamplePathFBM_Cholesky(double maturity, int numTimeSteps, double H)
{//Simulation of a fractional brownian motion sample path with exponent H by Cholesky decomposition
	double dt = maturity / numTimeSteps;
	default_random_engine generator(random_device{}());
	normal_distribution<double> distribution(0, 1);
	vector<double> gaussianSample(numTimeSteps, 0);
	cSquareMatrix covarianceMatrix(numTimeSteps);
	for (int i = 0; i < numTimeSteps; i++)
	{
		gaussianSample[i] = distribution(generator);
		for (int j = 0; j < numTimeSteps; j++)
			covarianceMatrix(i, j) = CovarianceFBM((i + 1.) * dt, (j + 1.) * dt, H);
	}
	unique_ptr<cSquareMatrix> sigma = covarianceMatrix.Cholesky();
	unique_ptr<vector<double>> res(*sigma * gaussianSample);
	res->insert(res->begin(), 0);
	return res;
}

double 
CosinusKL(double t, int k_index_KL, double maturityKL)
{//This function will be useful to simulate a fractional brownian motion sample path by 
 //Series Expansion method
	return cos(M_PI * t * k_index_KL / maturityKL);
}

double
SinusKL(double t, int k_index_KL, double maturityKL)
{//This function will be useful to simulate a fractional brownian motion sample path by 
 //Series Expansion method
	return sin(M_PI * t * k_index_KL / maturityKL);
}

double 
F_KL(double H, double maturityKL, int k_index_KL, double variableOfIntegration)
{//This function will be useful to simulate a fractional brownian motion sample path with Hurst exponent H
 //by Series Expansion method
	return pow(variableOfIntegration, 2 * H)* CosinusKL(variableOfIntegration, k_index_KL, maturityKL);
}

double
G_KL(double H, double maturityKL, int k_index_KL, double variableOfIntegration)
{//This function will be useful to simulate a fractional brownian motion sample path with Hurst exponent H
 //by Series Expansion method
	return pow(variableOfIntegration, 2 * H - 2) * CosinusKL(variableOfIntegration, k_index_KL, maturityKL);
}

double 
Coeff_C0_KL(double H, double maturityKL)
{//This function will be useful to simulate a fractional brownian motion sample path with Hurst exponent H
 //by Series Expansion method
	if (H < 0.5)
		return 0;
	return H * pow(maturityKL, 2 * H - 2);
}

double 
Coeff_Ck_KL(double H, double maturityKL, int k_index_KL)
{//This function will be useful to simulate a fractional brownian motion sample path with Hurst exponent H
 //by Series Expansion method
	if (H < 0.5)
		return (2 / maturityKL) 
		* IntegrateForSeriesExpansionMethod(F_KL, H, maturityKL, k_index_KL, 0, maturityKL);
	return (-4 * H * (2 * H - 1) * maturityKL / pow(k_index_KL * M_PI, 2))
		* IntegrateForSeriesExpansionMethod(G_KL, H, maturityKL, k_index_KL, 0, maturityKL);
}

unique_ptr<vector<double>> 
SamplePathFBM_KL(double maturity, int numTimeSteps, double H, int upperBoundInSumKL)
{//Simulation of a fractional brownian motion sample path with Hurst exponent H by Series Expansion method
	double dt = maturity / numTimeSteps;
	default_random_engine generator(random_device{}());
	normal_distribution<double> distribution(0, 1);
	vector<double> gaussianSample(2 * (size_t)upperBoundInSumKL + 1, 0);
	vector<double> coeff_Ck(upperBoundInSumKL, 0);
	unique_ptr<vector<double>> FBM(new vector<double>(numTimeSteps, 0));
	for (int k_indexKL = 0; k_indexKL < upperBoundInSumKL; k_indexKL++)
	{//We generate the gaussian sample and we compute the coefficients ck
		coeff_Ck[k_indexKL] = Coeff_Ck_KL(H, maturity, k_indexKL + 1);
		gaussianSample[2 * (size_t)k_indexKL] = distribution(generator);
		gaussianSample[2 * (size_t)k_indexKL + 1] = distribution(generator);
	}
	for (int timeStep = 0; timeStep < numTimeSteps; timeStep++) //We compute the sum for each time step
	{
		double res = sqrt(Coeff_C0_KL(H, maturity)) * (timeStep + 1.) * dt * gaussianSample[0];
		for (int k_indexKL = 0; k_indexKL < upperBoundInSumKL; k_indexKL++)
		{
			try //DEBUG
			{
				if (coeff_Ck[k_indexKL] >= 0) throw "sqrt(negative number)";
			}
			catch (const char* err)
			{
				cerr << err << endl;
				exit(1);
			}
			double temp = sqrt(-coeff_Ck[k_indexKL] / 2);
			res += temp * (SinusKL((timeStep + 1.) * dt, k_indexKL + 1, maturity) 
				* gaussianSample[2 * (size_t)k_indexKL + 1]
				+ (1 - CosinusKL((timeStep + 1.) * dt, k_indexKL + 1, maturity)) 
				* gaussianSample[2 * (size_t)k_indexKL + 2]);
		}
		(*FBM)[timeStep] = res;
	}
	FBM->insert(FBM->begin(), 0); //We insert 0 so that the trajectory starts at 0
	return FBM;
}

void 
ExportData(const unique_ptr<vector<double>>& data, string fileName)
{//Export a FBM sample path to a text file
	ofstream file(fileName, ios::out | ios::trunc);
	if (file)
	{
		for (unsigned int i = 0; i < data->size(); i++)
			file << (*data)[i] << ',';
		file.close();
	}
	else
	{
		cerr << "Error while opening file" << endl;
		exit(1);
	}
}

unique_ptr<vector<double>> 
StandardBM(int numTimeSteps, const unique_ptr<cSquareMatrix>& sigma)
{//Simulation of a standard brownian motion sample path with Cholesky decomposition matrix sigma as an argument
 //so we just have to compute it once when we use a Monte Carlo algorithm.
	default_random_engine generator(random_device{}());
	normal_distribution<double> distribution(0, 1);
	vector<double> gaussianSample(numTimeSteps, 0);
	for (int i = 0; i < numTimeSteps; i++)
		gaussianSample[i] = distribution(generator);
	unique_ptr<vector<double>> res(*sigma * gaussianSample);
	res->insert(res->begin(), 0); //We insert 0 so that the trajectory starts at 0
	return res;
}

double 
BlackScholesCallPrice(
	double maturity,
	double strike,
	double volatility,
	double shortRate,
	double S0)
{//Black Scholes formula for a call option
	double d1 = (log(S0 / strike) + (shortRate + volatility * volatility / 2) * maturity) /
		(volatility * sqrt(maturity));
	double d2 = d1 - volatility * sqrt(maturity);
	return S0 * (erfc(-d1 / sqrt(2)) / 2) - 
		strike * exp(-shortRate * maturity) * (erfc(-d2 / sqrt(2)) / 2);
}

double
DiffusionRoughHeston(
	int numTimeSteps, 
	double maturity, 
	double P0, //stock price at time 0
	double Y0, //variance at time 0
	double alpha, 
	double lambdaStar, 
	double mu, 
	double beta,
	double C,
	double phi1Norm,
	double phi2Norm,
	const unique_ptr<cSquareMatrix>& sigma)
{//Returns the approximate for the asset price P at maturity date in the rHeston model using an Euler scheme.
 //Y0 is the variance at time 0 and P0 is the price at time 0. sigma is the volatility matrix that we get by
 //Cholesky decomposition for the standard brownian motion.
	double dt = maturity / numTimeSteps;
	double Y = Y0; //current variance
	double P = P0; //current stock price
	double rho = (1 - beta) / sqrt(2 * (1 + beta * beta));
	double lambda = alpha * lambdaStar / (C * tgamma(1 - alpha));
	double C2 = lambda * sqrt((1 + beta * beta) / (lambdaStar * mu * (1 + beta)));
	double C3 = (1 / (1 - (phi1Norm - phi2Norm))) * sqrt(2 / (1 + beta));
	const unique_ptr<vector<double>>& B(StandardBM(numTimeSteps, sigma));
	const unique_ptr<vector<double>>& B_orthogonal(StandardBM(numTimeSteps, sigma));
	const unique_ptr<vector<double>>& W(*(rho * (*B)) + *(sqrt(1 - rho * rho) * (*B_orthogonal)));
	unique_ptr<vector<double>> Y_memory(new vector<double>);
	unique_ptr<vector<double>> P_memory(new vector<double>);
	Y_memory->push_back(Y);
	P_memory->push_back(P);
	size_t i;
	for (i = 0; i < numTimeSteps; i++)
	{
		double drift = 0;
		double vol = 0;
		for (size_t j = 0; j < Y_memory->size(); j++) //We approximate the integrals in the dynamics of Y
		{
			drift += (alpha - 1) * pow((i + 1.) * dt - j * dt, alpha - 2) * lambda
				* ((1 + beta) - (*Y_memory)[j]) * dt / tgamma(alpha);
			vol += (alpha - 1) * pow((i + 1.) * dt - j * dt, alpha - 2) * C2 
				* sqrt((*Y_memory)[j]) * ((*B)[j + 1] - (*B)[j]) / tgamma(alpha);
		}
		P += C3 * sqrt(Y) * ((*W)[i + 1] - (*W)[i]);
		Y = Y + drift * dt + vol * dt;
		Y_memory->push_back(Y);
		P_memory->push_back(P);
	}
	string name = "rHestonPrice.txt"; //We export the stock price trajectory
	ExportData(P_memory, name);
	return P;
}

double 
IntegrandBergomi(double x, double s, double gamma)
{//Integrand in the covariance function of Wtilde in the rBergomi model
	return 1 / pow((1 - s) * (x - s), gamma);
}

double 
IntegrateCovBergomi(
	double(*f)(double, double, double),
	double x,
	double gamma,
	double lowerBound,
	double upperBound)
{//Aproximation to an integral using rectangle method over [lowerBound, upperBound] used to 
 //simulate the rBergomi diffusion
	int numSteps = 10000;
	double dt = (upperBound - lowerBound) / numSteps;
	double res = 0;
	for (int i = 0; i < numSteps; i++)
	{
		res += dt * f(x, lowerBound + (i + 1/2) * dt, gamma);
	}
	return res;
}

double 
Covariance_Wtilde(double s, double t, double HurstExponent)
{//Covariance function of Wtilde in the rBergomi model
	double gamma = 0.5 - HurstExponent;
	double _max = max(s, t);
	double _min = min(s, t);
	return 2 * HurstExponent* pow(min(s, t), 2 * HurstExponent) 
		* IntegrateCovBergomi(IntegrandBergomi, _max / _min, gamma, 0, 1);
}

double 
Covariance_StandardBM(double s, double t)
{//Covariance of the standard brownian motion Z in the rBergomi model
	return min(s, t);
}

double 
Covariance_Wtilde_Z(double s, double t, double HurstExponent, double rho)
{//Covariance between Wtilde_s and Z_t in the rBergomi model
	double DH = sqrt(2 * HurstExponent) / (HurstExponent + 0.5);
	return rho * DH* (pow(s, HurstExponent + 0.5) - pow(s - min(s, t), HurstExponent + 0.5));
}

unique_ptr<vector<double>> 
Trajectory_Wtilde_and_Z(int numTimeSteps, const unique_ptr<cSquareMatrix>& sigma)
{//Sample path for Wtilde and Z in the rBergomi model. The first half of the vector corresponds to Z and the
 //second half corresponds to Wtilde. Both trajectories begin at 0, so we insert 0 at the beginning and at the
 //middle of the vector. sigma is the square root of the covariance matrix for Wtilde and Z.
	default_random_engine generator(random_device{}());
	normal_distribution<double> distribution(0, 1);
	vector<double> gaussianSample(2 * (size_t)numTimeSteps, 0);
	for (int i = 0; i < 2 * numTimeSteps; i++)
		gaussianSample[i] = distribution(generator);
	unique_ptr<vector<double>> res(*sigma * gaussianSample);
	res->insert(res->begin() + numTimeSteps, 0);
	res->insert(res->begin(), 0);
	return res;
}

double
DiffusionRoughBergomi(
	int numTimeSteps,
	double maturity,
	double S0, //stock price at time 0
	double v0, //variance at time 0
	double eta,
	const unique_ptr<cSquareMatrix>& sigma)
{//Returns the approximate for the asset price S at maturity date in the rBergomi model using an Euler scheme.
 //v0 is the variance at time 0 and S0 is the price at time 0. sigma is the volatility matrix that we get by
 //Cholesky decomposition for the correlated processes Wtilde and Z.
	double S = S0; //current stock price
	double v = v0; //current variance
	const unique_ptr<vector<double>>& 
		WtildeAndZ(Trajectory_Wtilde_and_Z(numTimeSteps, sigma)); //Size 2*numTimeSteps + 2
																  //The +2 corresponds to the value "0" we add 
									                              //at the beginning of both trajectories
	unique_ptr<vector<double>> S_memory(new vector<double>);
	S_memory->push_back(S);
	for (size_t i = 0; i < numTimeSteps; i++)
	{
		S *= (1 + sqrt(v) * ((*WtildeAndZ)[i + 1] - (*WtildeAndZ)[i]));
		v *= (1 + eta * ((*WtildeAndZ)[numTimeSteps + i + 2] - (*WtildeAndZ)[numTimeSteps + i + 1]));
		S_memory->push_back(S);
	}
	string name = "rBergomiPrice.txt"; //We export the stock price trajectory
	ExportData(S_memory, name);
	return S;
}

double 
MonteCarlo_rBergomi_Call(
	int numMonteCarloSimulations,
	int numTimeSteps,
	double maturity,
	double S0,
	double v0,
	double eta,
	const unique_ptr<cSquareMatrix>& sigma,
	double strike,
	double shortRate)
{
	double alpha = 5.0 / 100; //for confidence interval
	double res = 0;
	double variance = 0;
	vector<double> payoffSimulationTab(numMonteCarloSimulations, 0);
	for (int i = 0; i < numMonteCarloSimulations; i++)
	{
		double tmp = fmax(DiffusionRoughBergomi(numTimeSteps, maturity, S0, v0, eta, sigma) - strike, 0)
			* exp(-shortRate * maturity);
		res += tmp;
		payoffSimulationTab[i] = tmp;
	}
	res /= numMonteCarloSimulations;
	for (int i = 0; i < numMonteCarloSimulations; i++)
	{
		variance += pow(payoffSimulationTab[i] - res, 2);
	}
	variance /= (numMonteCarloSimulations - 1.0);
	//show confidence interval with alpha = 5%
	cout << "Confidence Interval : [" << res - sqrt(variance / numMonteCarloSimulations) * 1.96
		<< ", " << res + sqrt(variance / numMonteCarloSimulations) * 1.96 << "]" << endl;
	return res;
}

double 
Weight_LiftedHeston(int n, int i, double rn, double H)
{//Returns the weight coefficient c_i_n in the Lifted Heston model
	double alpha = H + 0.5;
	return (pow(rn, 1 - alpha) - 1) * pow(rn, (alpha - 1) * (1 + n / 2.0)) * pow(rn, (1 - alpha) * i) / 
		(tgamma(alpha) * tgamma(2 - alpha));
}

double 
ReversionSpeed_LiftedHeston(int n, int i, double rn, double H)
{//Returns the reversion speed coefficient x_i_n in the Lifted Heston model
	double alpha = H + 0.5;
	return (1 - alpha) * (pow(rn, 2 - alpha) - 1) * pow(rn, i - 1.0 - n / 2.0) / 
		((2 - alpha) * (pow(rn, 1 - alpha) - 1));
}

double 
InitialForwardVarianceCurve(double t, double V0, double theta, int n, double rn, double H)
{//Returns the value g0(t) of the initial forward variance curve at time t
 //V0 and theta are supposed to be calibrated using variance swaps from the market
	double sum = 0;
	for (int i = 1; i <= n; i++)
	{
		double integralApproximate = 0;
		double dt = t / 10;
		for (int j = 0; j < 10; j++)
		{
			integralApproximate += dt * exp(-ReversionSpeed_LiftedHeston(n, i, rn, H) * (t - j * dt));
		}
		sum += integralApproximate * Weight_LiftedHeston(n, i, rn, H);
	}
	sum *= theta;
	return V0 + sum;
}

double DiffusionLiftedHeston(
	int numTimeSteps,
	double maturity,
	int n, //number of factors
	double rn,
	double S0, //stock price at time 0
	double V0, //variance at time 0
	double theta,
	double rho,
	double nu,
	double lambda,
	double H,
	const unique_ptr<cSquareMatrix>& sigma) //volatility matrix for standard brownian motion
{//Returns the approximate for the asset price S at maturity date in the Lifted Heston model using the 
 //"explicit-implicit" scheme.
	double dt = maturity / numTimeSteps;
	double V = V0; //current variance
	double S = S0; //current stock price
	const unique_ptr<vector<double>> & W(StandardBM(numTimeSteps, sigma));
	const unique_ptr<vector<double>> & W_orthogonal(StandardBM(numTimeSteps, sigma));
	const unique_ptr<vector<double>> & B(*(rho * (*W)) + *(sqrt(1 - rho * rho) * (*W_orthogonal)));
	unique_ptr<vector<double>> U(new vector<double>(n, 0)); //U_1_t, ..., U_n_t
	unique_ptr<vector<double>> V_memory(new vector<double>);
	unique_ptr<vector<double>> S_memory(new vector<double>);
	V_memory->push_back(V);
	S_memory->push_back(S);
	for (size_t i = 0; i < numTimeSteps; i++)
	{
		for (size_t j = 1; j <= n; j++)
		{
			(*U)[j - 1] = ((*U)[j - 1] - lambda * fmax(V, 0) * dt + nu * sqrt(fmax(V, 0)) *
				((*W)[i + 1] - (*W)[i])) / (1 + ReversionSpeed_LiftedHeston(n, j, rn, H) * dt);
		}
		V = InitialForwardVarianceCurve(i * dt, V0, theta, n, rn, H);
		for (size_t j = 1; j <= n; j++)
		{
			V += Weight_LiftedHeston(n, j, rn, H) * (*U)[j - 1];
		}
		S *= (1 + sqrt(fmax(V, 0)) * ((*B)[i + 1] - (*B)[i]));
		//cout << V << endl;
		V_memory->push_back(V);
		S_memory->push_back(S);
	}
	string name = "LiftedHestonPrice.txt"; //We export the stock price trajectory
	ExportData(S_memory, name);
	return S;
}

unique_ptr<vector<double>> 
CalibrateLiftedHeston(
	const unique_ptr<vector<double>>& marketPrices,
	const unique_ptr<vector<double>>& maturities,
	const unique_ptr<vector<double>>& strikes,
	string financialProduct,
	int numMonteCarloSimulations,
	int numTimeSteps,
	const unique_ptr<cSquareMatrix>& sigma,
	double shortRate)
{//Approximate optimal parameters for the Lifted Heston model
	double min = 0;
	double rho_opti = -1;
	double H_opti = 0.05;
	double theta_opti = 0.01;
	double lambda_opti = 0.1;
	double nu_opti = 0.1;
	double V0_opti = 0.01;
	for (double rho = -1; rho <= 1; rho += 1)
		for (double H = 0.05; H <= 0.25; H += 0.1)
			for (double theta = 0.01; theta <= 0.05; theta += 0.02)
				for (double lambda = 0.1; lambda <= 0.5; lambda += 0.2)
					for (double nu = 0.1; nu <= 0.5; nu += 0.2)
						for (double V0 = 0.01; V0 <= 0.05; V0 += 0.02)
						{
							double sum = 0;
							for (int i = 0; i < maturities->size(); i++)
								for (int j = 0; j < strikes->size(); j++)
								{
									int k = i * strikes->size() + j;
									sum += pow((*marketPrices)[k] - MonteCarlo_pricing(
										financialProduct, "LiftedHeston", numMonteCarloSimulations,
										numTimeSteps, (*maturities)[i], sigma, (*strikes)[j],
										shortRate), 2);
								}
							if (sum < min || min == 0)
							{
								min = sum;
								rho_opti = rho;
								theta_opti = theta;
								H_opti = H;
								lambda_opti = lambda;
								nu_opti = nu;
								V0_opti = V0;
							}
						}
	return unique_ptr<vector<double>>(new vector<double>
		{ V0_opti,theta_opti,lambda_opti,nu_opti,rho_opti,H_opti });
}

double 
ImplicitVolLiftedHeston(
	double LiftedHestonPrice,
	double maturity,
	double strike,
	double shortRate)
{//Returns an approximate for the volatility in the BS model which matches the Lifted Heston call price
	double x = 0;
	double y = 1;
	for (int i = 0; i < 100; i++)
	{
		double z = (x + y) / 2;
		double BlackScholesPrice = BlackScholesCallPrice(maturity, strike, z / (1 - z), shortRate, 2930);
		if (LiftedHestonPrice > BlackScholesPrice)
			x = z;
		else
			y = z;
	}
	//cout << LiftedHestonPrice << endl;
	double z = (x + y) / 2;
	//cout << BlackScholesCallPrice(maturity, strike, z / (1 - z), shortRate, 2930) << endl;
	//cout << z / (1 - z) << endl;
	return z / (1 - z);
}

unique_ptr<vector<double>> 
GenerateVolSurf_LiftedHeston(
	double shortRate,
	const unique_ptr<cSquareMatrix>& sigma)
{//Generate volatility surface from the Lifted Heston model
 //The result contains a vector as : sigma(T1,K1), ... ,sigma(T1,Kn),sigma(T2,K1), ...... ,sigma(Tn,Kn)
 //There are 24 maturities ans 19 strikes
	unique_ptr<vector<double>> volSurface(new vector<double>(19 * 24, 0));
	size_t i;
	size_t j;
	vector<double> maturity{ 1.0 / 12, 2.0 / 12, 3.0 / 12, 4.0 / 12, 5.0 / 12, 6.0 / 12,
		7.0 / 12, 8.0 / 12, 9.0 / 12, 10.0 / 12, 11.0 / 12, 12.0 / 12, 13.0 / 12, 14.0 / 12,
		15.0 / 12, 16.0 / 12, 17.0 / 12, 18.0 / 12, 19.0 / 12, 20.0 / 12, 21.0 / 12,
		22.0 / 12, 23.0 / 12, 24.0 / 12 };
	vector<double> strike{ 5.0 * 2930 / 10, 6.0 * 2930 / 10, 7.0 * 2930 / 10, 7.5 * 2930 / 10,
		8.0 * 2930 / 10, 8.5 * 2930 / 10, 9.0 * 2930 / 10, 9.5 * 2930 / 10, 9.75 * 2930 / 10,
		10.0 * 2930 / 10, 10.25 * 2930 / 10, 10.5 * 2930 / 10, 11.0 * 2930 / 10,
		11.5 * 2930 / 10, 12.0 * 2930 / 10, 12.5 * 2930 / 10, 13.0 * 2930 / 10, 14.0 * 2930 / 10,
		15.0 * 2930 / 10 };
	for (i = 0; i < 24; i++)
	{
		for (j = 0; j < 19; j++)
		{
			double LiftedHestonPrice = MonteCarlo_pricing("Call", "LiftedHeston",
				120, 15, maturity[i], sigma, strike[j], shortRate);
			//cout << "(" << maturity[i] << ", " << strike[j] << ") " <<
				//"LiftedPrice : " << LiftedHestonPrice << endl;
			(*volSurface)[i * 19 + j] = ImplicitVolLiftedHeston(LiftedHestonPrice,
				maturity[i], strike[j], shortRate);
			//cout << (*volSurface)[i * 19 + j] << endl;
		}
	}
	string name = "LiftedHestonVolSurface.txt";
	ExportData(volSurface, name);
	return volSurface;
}

double MonteCarlo_pricing(
	string financialProduct,
	string model,
	int numMonteCarloSimulations,
	int numTimeSteps,
	double maturity,
	const unique_ptr<cSquareMatrix>& sigma,
	double strike,
	double shortRate)
{
	double alpha = 5.0 / 100; //for confidence interval
	double res = 0;
	double variance = 0;
	vector<double> payoffSimulationTab(numMonteCarloSimulations, 0);
	for (int i = 0; i < numMonteCarloSimulations; i++)
	{
		double tmp = 0;
		if (financialProduct == "Call" && model == "rBergomi")
			tmp = fmax(DiffusionRoughBergomi(numTimeSteps, maturity, 100, 0.005, 0.2, sigma) - strike, 0) 
			* exp(-shortRate * maturity);
		if (financialProduct == "Put" && model == "rBergomi")
			tmp = fmax(strike - DiffusionRoughBergomi(numTimeSteps, maturity, 100, 0.005, 0.2, sigma), 0)
			* exp(-shortRate * maturity);
		if (financialProduct == "Call" && model == "LiftedHeston")
			tmp = fmax(DiffusionLiftedHeston(numTimeSteps, maturity, 20, 2.5, 2930, 0.02, 0.2,
				-0.7, 0.3, 0.3, 0.1, sigma) - strike, 0) * exp(-shortRate * maturity);
		if (financialProduct == "Put" && model == "LiftedHeston")
			tmp = fmax(strike - DiffusionLiftedHeston(numTimeSteps, maturity, 20, 2.5, 100, 0.02, 0.02,
				-0.7, 0.3, 0.3, 0.1, sigma), 0) * exp(-shortRate * maturity);
		res += tmp;
		payoffSimulationTab[i] = tmp;
	}
	res /= numMonteCarloSimulations;
	/*for (int i = 0; i < numMonteCarloSimulations; i++)
	{
		variance += pow(payoffSimulationTab[i] - res, 2);
	}
	variance /= (numMonteCarloSimulations - 1.0);*/
	//show confidence interval with alpha = 5%
	//cout << "Confidence Interval : [" << res - sqrt(variance / numMonteCarloSimulations) * 1.96
		//<< ", " << res + sqrt(variance / numMonteCarloSimulations) * 1.96 << "]" << endl;
	//cout << res << endl;
	return res;
}

