#include "ThomasAlgorithm.h"

#include <vector>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

using std::cout;
using std::endl;
using boost::numeric::ublas::matrix;
using boost::numeric::ublas::identity_matrix;
using boost::numeric::ublas::zero_matrix;
using boost::numeric::ublas::scalar_matrix;
using boost::numeric::ublas::prod;

void testThomasAlgorithm()
{
	std::vector<double> A =
	{
		3, -1, 0, 0, 0,
		1, -4, 2, 0, 0,
		0, -1, 4, 2, 0,
		0, 0, 3, 5, -1,
		0, 0, 0, 1, 2
	};

	std::vector<double> d =
	{
		1, 2, 3, 4, 5,
	};

	ThomasAlgorithm ta(A, d);
	ThomasAlgorithm::Vector x = ta.solve();
	cout << "solution vector x is " << endl << x << endl;
}

ThomasAlgorithm::ThomasAlgorithm()
	: A(Matrix(0,0)), d(Vector(0,0))
{
}

ThomasAlgorithm::ThomasAlgorithm(const Matrix& A, const Vector& d)
	: A(A), d(d)
{
}

ThomasAlgorithm::ThomasAlgorithm(
	const std::vector<double>& A, const std::vector<double>& d)
{
	int eqtSize = d.size();

	Matrix A_(eqtSize, eqtSize);
	for (int i = 0; i < eqtSize; ++i)
		for (int j = 0; j < eqtSize; ++j)
			A_(i, j) = A[i * eqtSize + j];
	this->A = A_;

	Vector d_(eqtSize);
	for (int i = 0; i < eqtSize; ++i)
		d_(i) = d[i];
	this->d = d_;
}

ThomasAlgorithm::ThomasAlgorithm(const ThomasAlgorithm& ta)
	: A(ta.A), d(ta.d)
{
}

ThomasAlgorithm& ThomasAlgorithm::operator=(const ThomasAlgorithm& ta)
{
	if (this != &ta) {
		A = ta.A;
		d = ta.d;
	}
	return *this;
}

ThomasAlgorithm::~ThomasAlgorithm()
{
}

int ThomasAlgorithm::eqtSize()
{
	return A.size1();
}

// Core Function
ThomasAlgorithm::Vector ThomasAlgorithm::solve()
{
	int size = eqtSize();
	Vector a_(size);	// a'
	Vector d_(size);	// d'
	Vector x(size); // x (solution)

	// calculate a' and d'
	a_(0) = A(0, 0);
	d_(0) = d(0);
	for (int i = 1; i < size; ++i) {
		a_(i) = A(i, i) - A(i, i - 1) * A(i - 1, i) / a_(i - 1);
		d_(i) = d(i) - A(i, i - 1) * d_(i - 1) / a_(i - 1);
	}

	// calculate x
	x(size - 1) = d_(size - 1) / a_(size - 1);
	for (int i = size - 2; i >= 0; --i) {
		x(i) = d_(i) / a_(i) - A(i, i + 1) / a_(i) * x(i + 1);
	}

	return x;
}

bool ThomasAlgorithm::isSquare()
{
	// 미완성
	return true;
}

bool ThomasAlgorithm::isTridiagonal()
{
	// 미완성
	return true;
}

bool ThomasAlgorithm::isStrictlyDiagonallyDominant()
{
	// 미완성
	return true;
}

