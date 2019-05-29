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
	
	int n = 100;
	double maturity = 1;
	double dt = maturity / n;
	cSquareMatrix CovarianceStandardBM(n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			CovarianceStandardBM(i, j) = CovarianceFBM((i + 1.) * dt, (j + 1.) * dt, 0.5);
	unique_ptr<cSquareMatrix> sigma = CovarianceStandardBM.Cholesky();
	for (int i = 0; i < 100; i++)
	{
		double stockPrice = DiffusionRoughHeston(n, maturity, 100, 5, 0.65, 1, 1, 1.2, 1, 0.5, 0.25, sigma);
		cout << stockPrice << endl;
	}
	//cout << stockPrice << endl;
	

}