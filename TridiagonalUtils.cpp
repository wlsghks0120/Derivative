#include "TridiagonalUtils.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

typedef boost::numeric::ublas::matrix<double> Matrix;
typedef boost::numeric::ublas::vector<double> Vector;

using std::cout;
using std::endl;

Vector initVector(const std::vector<double>& vec)
{
	Vector res(vec.size());
	for (int i = 0; i < (int)vec.size(); ++i)
		res(i) = vec[i];
	return res;
}

Matrix createMatrix(const std::vector<std::vector<double>>& mat)
{
	Matrix res(mat.size(), mat[0].size());
	for (int i = 0; i < (int)mat.size(); ++i)
		for (int j = 0; j < (int)mat[0].size(); ++j)
			res(i, j) = mat[i][j];
	return res;
}

// a : diagonal vector
// b,c : subdiagonal vector
Matrix createTridiagonalMatrix(const Vector& a, const Vector& b, const Vector& c)
{
	int size = a.size();
	Matrix res(size, size);
	for (int i = 0; i < size; ++i)
		for (int j = 0; j < size; ++j)
			res(i, j) = 0;
	for (int i = 0; i < size; ++i)
		res(i, i) = a(i);
	for (int i = 0; i < size - 1; ++i)
		res(i + 1, i) = b(i);
	for (int i = 0; i < size - 1; ++i)
		res(i, i + 1) = c(i);
	return res;
}

void printMatrix(const Matrix& mat)
{
	cout << endl;
	for (int i = 0; i < mat.size1(); ++i) {
		for (int j = 0; j < mat.size2(); ++j) {
			cout << mat(i, j) << ", ";
		}
		cout << endl;
	}
}

void printVector(const Vector& vec)
{
	cout << endl;
	for (int i = 0; i < vec.size(); ++i)
		cout << vec(i) << ",";
	cout << endl;
}