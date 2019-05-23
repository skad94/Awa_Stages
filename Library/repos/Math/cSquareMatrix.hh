#ifndef CSQUAREMATRIX_HH
#define CSQUAREMATRIX_HH

#include <vector>
#include <cmath>
#include <math.h>
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
	unique_ptr<cSquareMatrix> operator*(const cSquareMatrix& m) const;
	unique_ptr<vector<double>> operator*(const vector<double>& v) const;
	unique_ptr<cSquareMatrix> Transpose() const;
	unique_ptr<cSquareMatrix> Cholesky() const;
};

void ShowVector(const vector<double>& v);

#endif