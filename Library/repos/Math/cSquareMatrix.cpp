#include "cSquareMatrix.hh"


cSquareMatrix::cSquareMatrix(const int& size)
	:_size(size), _data(_size* _size, 0) {}

cSquareMatrix::cSquareMatrix(const cSquareMatrix& m)
	: _size(m._size), _data(m._data) {}

cSquareMatrix& 
cSquareMatrix::operator=(const cSquareMatrix& m)
{
	_size = m._size;
	_data = m._data;
	return *this;
}

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
cSquareMatrix::ShowMatrix() const
{//Display contents of a square matrix
	for (int i = 0; i < _size; i++)
	{
		for (int j = 0; j < _size; j++)
		{
			cout << operator()(i, j) << "   ";
		}
		cout << endl;
	}
}

cSquareMatrix& 
cSquareMatrix::operator*=(const cSquareMatrix& m)
{//Matrix multiplication
	cSquareMatrix res(*this);
	double s;
	for (int i = 0; i < _size; i++)
	{
		for (int j = 0; j < _size; j++)
		{
			s = 0;
			for (int k = 0; k < _size; k++)
			{
				s += res(i, k) * m(k, j);
			}
			operator()(i, j) = s;
		}
	}
	return *this;
}

cSquareMatrix 
cSquareMatrix::operator*(const cSquareMatrix& m) const
{
	cSquareMatrix res(*this);
	res *= m;
	return res;
}

vector<double> 
cSquareMatrix::operator*(const vector<double>& v) const
{//Product of matrix and vector
	vector<double> res(_size, 0);
	double s;
	for (int i = 0; i < _size; i++)
	{
		s = 0;
		for (int k = 0; k < _size; k++)
		{
			s += operator()(i, k) * v[k];
		}
		res[i] = s;
	}
	return res;
}

cSquareMatrix 
cSquareMatrix::Transpose()
{//Matrix transposition
	cSquareMatrix res(_size);
	for (int i = 0; i < _size; i++)
	{
		for (int j = 0; j < _size; j++)
		{
			res(i, j) = operator()(j, i);
		}
	}
	return res;
}

cSquareMatrix 
cSquareMatrix::Cholesky()
{//Cholesky decomposition
	cSquareMatrix L(_size);
	double s;
	for (int i = 0; i < _size; i++)
	{
		for (int j = i; j < _size; j++)
		{
			s = 0;
			for (int k = 0; k < i; k++)
			{
				s += L(i, k) * L(j, k);
			}
			if (i == j) L(i, j) = sqrt(operator()(i, j) - s);
			else L(j, i) = (operator()(i, j) - s) / L(i, i);
		}
	}
	return L;
}

void 
ShowVector(const vector<double>& v) 
{//Display contents of a column vector
    for (int i = 0; i < v.size(); i++) {
        cout << v[i] << endl;
    }
}


