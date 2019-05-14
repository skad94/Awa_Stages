#include <vector>
#include <cmath>
#include <iostream>
using namespace std;

class cSquareMatrix {
    private:
    int _size;
    vector<double> _data;

    public:
    cSquareMatrix(int size);
    cSquareMatrix(const cSquareMatrix& m);
    ~cSquareMatrix();
    double& operator()(int i, int j);
    double operator()(int i, int j) const;
    int GetSize() const;
};

void ShowMatrix(const cSquareMatrix& m);
void ShowVector(const vector<double>& v);
cSquareMatrix operator*(const cSquareMatrix& m1, const cSquareMatrix& m2);
vector<double> operator*(const cSquareMatrix& m, const vector<double>& v);
cSquareMatrix Transpose(const cSquareMatrix& m);
cSquareMatrix Cholesky(const cSquareMatrix& m);