#pragma once

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <vector>

typedef boost::numeric::ublas::matrix<double> Matrix;
typedef boost::numeric::ublas::vector<double> Vector;

Vector createVector(const std::vector<double>& vec);
Matrix createMatrix(const std::vector<std::vector<double>>& mat);
Matrix createTridiagonalMatrix(const Vector& a, const Vector& b, const Vector& c);
void printMatrix(const Matrix& mat);
void printVector(const Vector& vec);