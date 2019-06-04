#define _USE_MATH_DEFINES 
#include "Functions.hh"

double 
IntegrateForSeriesExpansionMethod(double(*f)(double, double, int, double), 
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
	vector<double> gaussianSample(2 * upperBoundInSumKL + 1, 0);
	vector<double> coeff_Ck(upperBoundInSumKL, 0);
	unique_ptr<vector<double>> FBM(new vector<double>(numTimeSteps, 0));
	for (int k_indexKL = 0; k_indexKL < upperBoundInSumKL; k_indexKL++)
	{//We generate the gaussian sample and we compute the coefficients ck
		coeff_Ck[k_indexKL] = Coeff_Ck_KL(H, maturity, k_indexKL + 1);
		gaussianSample[2 * k_indexKL] = distribution(generator);
		gaussianSample[2 * k_indexKL + 1] = distribution(generator);
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
				* gaussianSample[2 * k_indexKL + 1]
				+ (1 - CosinusKL((timeStep + 1.) * dt, k_indexKL + 1, maturity)) 
				* gaussianSample[2 * k_indexKL + 2]);
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
DiffusionRoughHeston(
	int n, 
	double maturity, 
	double P0, 
	double Y0, 
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
	double dt = maturity / n;
	double Y = Y0;
	double P = P0;
	double rho = (1 - beta) / sqrt(2 * (1 + beta * beta));
	double lambda = alpha * lambdaStar / (C * tgamma(1 - alpha));
	double C2 = lambda * sqrt((1 + beta * beta) / (lambdaStar * mu * (1 + beta)));
	double C3 = (1 / (1 - (phi1Norm - phi2Norm))) * sqrt(2 / (1 + beta));
	const unique_ptr<vector<double>>& B(StandardBM(n, sigma));
	const unique_ptr<vector<double>>& B_orthogonal(StandardBM(n, sigma));
	const unique_ptr<vector<double>>& W(*(rho * (*B)) + *(sqrt(1 - rho * rho) * (*B_orthogonal)));
	unique_ptr<vector<double>> Y_memory(new vector<double>);
	unique_ptr<vector<double>> P_memory(new vector<double>);
	Y_memory->push_back(Y);
	P_memory->push_back(P);
	int i;
	for (i = 0; i < n; i++)
	{
		double drift = 0;
		double vol = 0;
		for (unsigned int j = 0; j < Y_memory->size(); j++)
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
	string name = "rHestonPrice.txt";
	ExportData(P_memory, name);
	return P;
}

double 
IntegrandBergomi(double x, double s, double gamma)
{//Integrand in the covariance function of Wtilde in the rBergomi model
	return 1 / pow((1 - s) * (x - s), gamma);
}

double 
IntegrateCovBergomi(double(*f)(double, double, double), double x, double gamma, double a, double b)
{//Aproximation to an integral using rectangle method over [a, b], used to simulate the rBergomi diffusion
	int n = 10000;
	double dt = (b - a) / n;
	double res = 0;
	for (int i = 0; i < n; i++)
	{
		res += dt * f(x, a + (i + 1/2) * dt, gamma);
	}
	return res;
}

double 
Covariance_Wtilde(double u, double v, double H)
{//Covariance function of Wtilde in the rBergomi model
	double gamma = 0.5 - H;
	double _max = max(u, v);
	double _min = min(u, v);
	return 2 * H * pow(min(u, v), 2 * H) * IntegrateCovBergomi(IntegrandBergomi, _max / _min, gamma, 0, 1);
}

double 
Covariance_StandardBM(double u, double v)
{//Covariance of the standard brownian motion Z in the rBergomi model
	return min(u, v);
}

double 
Covariance_Wtilde_Z(double u, double v, double H, double rho)
{//Covariance between Wtilde_u and Z_v in the rBergomi model
	double DH = sqrt(2 * H) / (H + 0.5);
	return rho * DH* (pow(u, H + 0.5) - pow(u - min(u, v), H + 0.5));
}

unique_ptr<vector<double>> 
Trajectory_Wtilde_and_Z(int n, const unique_ptr<cSquareMatrix>& sigma)
{//Sample path for Wtilde and Z in the rBergomi model. The first half of the vector corresponds to Z and the
 //second half corresponds to Wtilde. Both trajectories begin at 0, so we insert 0 at the beginning and at the
 //middle of the vector.
	default_random_engine generator(random_device{}());
	normal_distribution<double> distribution(0, 1);
	vector<double> gaussianSample(2 * n, 0);
	for (int i = 0; i < 2 * n; i++)
		gaussianSample[i] = distribution(generator);
	unique_ptr<vector<double>> res(*sigma * gaussianSample);
	res->insert(res->begin() + n, 0);
	res->insert(res->begin(), 0);
	return res;
}

double 
DiffusionRoughBergomi(
	int n,
	double maturity,
	double S0,
	double v0,
	double eta,
	const unique_ptr<cSquareMatrix>& sigma)
{//Returns the approximate for the asset price S at maturity date in the rBergomi model using an Euler scheme.
 //v0 is the variance at time 0 and S0 is the price at time 0. sigma is the volatility matrix that we get by
 //Cholesky decomposition for the correlated processes Wtilde and Z.
	double S = S0;
	double v = v0;
	const unique_ptr<vector<double>>& WtildeAndZ(Trajectory_Wtilde_and_Z(n, sigma)); //size 2*n + 2
	unique_ptr<vector<double>> S_memory(new vector<double>);
	S_memory->push_back(S);
	for (int i = 0; i < n; i++)
	{
		S *= (1 + sqrt(v) * ((*WtildeAndZ)[i + 1] - (*WtildeAndZ)[i]));
		v *= (1 + eta * ((*WtildeAndZ)[n + i + 2] - (*WtildeAndZ)[n + i + 1]));
		S_memory->push_back(S);
	}
	string name = "rBergomiPrice.txt";
	ExportData(S_memory, name);
	return S;
}