#pragma once

#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

/*
���� ������ġ�� �ȸ��������. std::vector �޴� �����ڿ���
A�� d�� ����� ��ġ�ϴ��� Ȯ���� ���ϰ� �ְ�
����� Tridiagonal ����, Square����, Strictly diagonally dominant ���� Ȯ�� �ؾ���.
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
