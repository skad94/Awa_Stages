#include <vector>
#include <cmath>
#include <iostream>
using namespace std;

class cSquareMatrix 
{
    private:
    int _size;
    vector<double> _data;

    public:
    cSquareMatrix(const int& size);
    cSquareMatrix(const cSquareMatrix& m);
	cSquareMatrix& operator=(const cSquareMatrix& m);
    double& operator()(const int& i, const int& j);
    double operator()(const int& i, const int& j) const;
    int GetSize() const;
	void ShowMatrix() const;
	cSquareMatrix& operator*=(const cSquareMatrix& m);
	cSquareMatrix operator*(const cSquareMatrix& m) const;
	vector<double> operator*(const vector<double>& v) const;
	cSquareMatrix Transpose();
	cSquareMatrix Cholesky();
};

void ShowVector(const vector<double>& v);