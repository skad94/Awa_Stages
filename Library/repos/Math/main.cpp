#include "cSquareMatrix.hh"


int main() {
    cSquareMatrix m1(3);
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
	m2.ShowMatrix();
}