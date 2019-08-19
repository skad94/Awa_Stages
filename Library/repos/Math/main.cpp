#include "Functions.hh"

int main() 
{
	srand((unsigned int)time(NULL));

    /*cSquareMatrix m1(3);
    m1(0,0) = 4;
    m1(0,1) = 2;
    m1(0,2) = 1;
    m1(1,0) = 2;
    m1(1,1) = 4;
    m1(1,2) = 1;
    m1(2,0) = 1;
    m1(2,1) = 1;
    m1(2,2) = 4;
    cSquareMatrix m2(3);
    m2(0,0) = 0;
    m2(0,1) = 1;
    m2(0,2) = 0;
    m2(1,0) = 1;
    m2(1,1) = 0;
    m2(1,2) = 2;
    m2(2,0) = 0;
    m2(2,1) = 0;
    m2(2,2) = 1;
    vector<double> v(3,0);
    v[0] = 1;
    v[1] = 0;
    v[2] = -1;
    (m2*m1).ShowMatrix();
	cout << endl;
	cSquareMatrix L = m1.Cholesky();
	L.ShowMatrix();
	cout << endl;
    cSquareMatrix LT = L.Transpose();
	LT.ShowMatrix();
	cout << endl;
    (L*LT).ShowMatrix(); // == m1
	cout << endl;
	m2 = L;
	m2.ShowMatrix();*/
	
	//unique_ptr<vector<double>> v(SamplePathFBM_Cholesky(3, 100, 0.3));
	//ShowVector(*v);

	/*const unique_ptr<vector<double>>& v1(SamplePathFBM_KL(3, 200, 0.8, 200));
	const unique_ptr<vector<double>>& v2(SamplePathFBM_Cholesky(3, 200, 0.8));
	const unique_ptr<vector<double>>& v3(SamplePathFBM_KL(3, 200, 0.3, 200));
	const unique_ptr<vector<double>>& v4(SamplePathFBM_Cholesky(3, 200, 0.3));
	string name1 = "SimulationFBM_KL_08.txt";
	string name2 = "SimulationFBM_Cholesky_08.txt";
	string name3 = "SimulationFBM_KL_03.txt";
	string name4 = "SimulationFBM_Cholesky_03.txt";
	ExportData(v1, name1);
	ExportData(v2, name2);
	ExportData(v3, name3);
	ExportData(v4, name4);*/



	//unique_ptr<cSquareMatrix> ;
	


	int numTimeSteps = 10;
	double maturity = 1;
	double dt = maturity / numTimeSteps;
	double HurstExponent = 0.1;
	double rho = -0.5;

	//rHeston
	/*cSquareMatrix CovarianceStandardBM(numTimeSteps); //We fill the covariance matrix for the standard BM
	for (int i = 0; i < numTimeSteps; i++)
		for (int j = 0; j < numTimeSteps; j++)
			CovarianceStandardBM(i, j) = Covariance_StandardBM((i + 1.) * dt, (j + 1.) * dt);
	unique_ptr<cSquareMatrix> sigmaStandard = CovarianceStandardBM.Cholesky();
	double stockPriceHeston = DiffusionRoughHeston(numTimeSteps, maturity,
		100, 10, 0.55, 1, 1, 1.2, 1, 0.5, 0.25, sigmaStandard);*/

	//rBergomi
	/*cSquareMatrix CovarianceMatrixBergomi(2 * numTimeSteps); //We fill the covariance matrix for Wtilde and Z
	for (int i = 0; i < 2 * numTimeSteps; i++)
		for (int j = 0; j < 2 * numTimeSteps; j++)
			if (i < numTimeSteps && j < numTimeSteps)
				CovarianceMatrixBergomi(i, j) = Covariance_StandardBM((i + 1.) * dt, (j + 1.) * dt);
			else if (i >= numTimeSteps && j >= numTimeSteps)
				CovarianceMatrixBergomi(i, j) = Covariance_Wtilde((i + 1.) * dt, (j + 1.) * dt, HurstExponent);
			else
				CovarianceMatrixBergomi(i, j) = Covariance_Wtilde_Z((i + 1.) * dt, (j + 1.) * dt,
					HurstExponent, rho);
	unique_ptr<cSquareMatrix> sigmaBergomi = CovarianceMatrixBergomi.Cholesky();*/
	//double stockPriceBergomi = DiffusionRoughBergomi(numTimeSteps, maturity, 100, 0.005, 0.2, sigmaBergomi);
	
	//Monte Carlo rBergomi
	/*double MCprice = MonteCarlo_rBergomi_Call(1000, numTimeSteps, maturity, 100, 0.005, 0.2,
		sigmaBergomi, 95, 0.01);
	cout << "Monte Carlo price : " << MCprice << endl;*/

	//Lifted Heston
	cSquareMatrix CovarianceStandardBM(numTimeSteps); //We fill the covariance matrix for the standard BM
	for (int i = 0; i < numTimeSteps; i++)
		for (int j = 0; j < numTimeSteps; j++)
			CovarianceStandardBM(i, j) = Covariance_StandardBM((i + 1.) * dt, (j + 1.) * dt);
	unique_ptr<cSquareMatrix> sigmaStandard = CovarianceStandardBM.Cholesky();
	/*double stockPriceLiftedHeston = DiffusionLiftedHeston(numTimeSteps, maturity, 20, 2.5, 10, 0.02,
		0.02, -0.7, 0.3, 0.3, HurstExponent, sigmaStandard);
	cout << stockPriceLiftedHeston << endl;
	unique_ptr<vector<double>> marketPrices(new vector<double>(1, 6));
	unique_ptr<vector<double>> maturities(new vector<double>(1, 1));
	unique_ptr<vector<double>> strikes(new vector<double>(1, 95));
	string financialProduct = "Call";
	ShowVector(*CalibrateLiftedHeston(marketPrices,
		maturities, strikes,
		financialProduct, 10, numTimeSteps,
		sigmaStandard, 0.01));*/

	/*unique_ptr<vector<double>> volatilitySurface = GenerateVolSurf_LiftedHeston(
		0, sigmaStandard);
	ShowVector(*volatilitySurface);*/
	
	vector<int> nn{ 5, 20, 100, 1000, 10000 };
	cout << MonteCarlo_pricing("Call", "rHeston", 500, 10, 1, sigmaStandard, 2900, 0) << endl;
	for (int i = 0; i < 5; i++)
		cout << MonteCarlo_pricing("Call", "LiftedHeston", 500, 10, 1, sigmaStandard, 2900, 0,
			nn[i], 1.0 + 10 / pow(nn[i], 0.9)) << endl;




}