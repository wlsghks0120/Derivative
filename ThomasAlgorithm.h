#pragma once

#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

/*
아직 안전장치가 안만들어졌음. std::vector 받는 생성자에서
A와 d의 사이즈가 일치하는지 확인을 안하고 있고
행렬이 Tridiagonal 인지, Square인지, Strictly diagonally dominant 인지 확인 해야함.
*/

class ThomasAlgorithm {
public:
	typedef boost::numeric::ublas::matrix<double> Matrix;
	typedef boost::numeric::ublas::vector<double> Vector;
	
	ThomasAlgorithm();
	ThomasAlgorithm(const Matrix& A, const Vector& d);
	ThomasAlgorithm(const std::vector<double>& A, const std::vector<double>& d);
	ThomasAlgorithm(const ThomasAlgorithm& ta);
	ThomasAlgorithm& operator=(const ThomasAlgorithm& ta);
	~ThomasAlgorithm();

	int eqtSize();
	Vector solve();

private:
	Matrix A;
	Vector d;

	bool isSquare();
	bool isTridiagonal();
	bool isStrictlyDiagonallyDominant();
};

void testThomasAlgorithm();
