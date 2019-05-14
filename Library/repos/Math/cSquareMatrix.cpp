#include "cSquareMatrix.hh"

cSquareMatrix::cSquareMatrix(const int& size)
	:_size(size), _data(_size* _size, 0) {}

cSquareMatrix::cSquareMatrix(const cSquareMatrix& m)
	: _size(m._size), _data(m._data) {}


double& 
cSquareMatrix::operator()(const int& i, const int& j) 
{
    return _data[i*_size + j];
}

double 
cSquareMatrix::operator()(const int& i, const int& j) const 
{
    return _data[i*_size + j];
}

int 
cSquareMatrix::GetSize() const 
{
    return _size;
}

void 
ShowMatrix(const cSquareMatrix& m) 
{//Display contents of a square matrix
    int n = m.GetSize();
    for (int i = 0; i < n; i++) 
	{
        for (int j = 0; j < n; j++) 
		{
            cout << m(i,j) << "   ";
        }
        cout << endl;
    }
}

void 
ShowVector(const vector<double>& v) 
{//Display contents of a column vector
    for (int i = 0; i < v.size(); i++) {
        cout << v[i] << endl;
    }
}

cSquareMatrix 
operator*(const cSquareMatrix& m1, const cSquareMatrix& m2) 
{//Matrix multiplication
    int n = m1.GetSize();
    cSquareMatrix m(n);
    double s;
    for (int i = 0; i < n; i++) 
	{
        for (int j = 0; j < n; j++) 
		{
            s = 0;
            for (int k = 0; k < n; k++) 
			{
                s += m1(i,k)*m2(k,j);
            }
            m(i,j) = s;
        }
    }
    return m;
}

vector<double> 
operator*(const cSquareMatrix& m, const vector<double>& v) 
{//Product of matrix and vector
    int n = v.size();
    vector<double> res(n,0);
    double s;
    for (int i = 0; i < n; i++) 
	{
        s = 0;
        for (int k = 0; k < n; k++) 
		{
            s += m(i,k)*v[k];
        }
        res[i] = s;
    }
    return res;
}

cSquareMatrix 
Transpose(const cSquareMatrix& m)
{//Matrix transposition
    int n = m.GetSize();
    cSquareMatrix res(n);
    for (int i = 0; i < n; i++) 
	{
        for (int j = 0; j < n; j++) 
		{
            res(i,j) = m(j,i);
        }
    }
    return res;
}

cSquareMatrix
Cholesky(const cSquareMatrix& m) 
{//Cholesky decomposition, return L such that "L LT = m"
    int n = m.GetSize();
    cSquareMatrix L(n);
    double s;
    for (int i = 0; i < n; i++) 
	{
        for (int j = i; j < n; j++) 
		{
            s = 0;
            for (int k = 0; k < i; k++) 
			{
                s += L(i,k)*L(j,k);
            }
            if (i == j) L(i,j) = sqrt(m(i,j) - s);
            else L(j,i) = (m(i,j) - s)/L(i,i);
        }
    }
    return L;
}