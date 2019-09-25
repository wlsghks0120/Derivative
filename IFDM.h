#pragma once

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

// �ϴ� european put option���� �ۼ�.
// ���߿� �ٸ� �Ļ���ǰ�鵵 Ȱ���� �� �ֵ��� Ȯ���غ���.

class IFDM {
public:
	typedef boost::numeric::ublas::vector<double> Vector;
	typedef boost::numeric::ublas::matrix<double> Matrix;

	IFDM();
	IFDM(int K, double r, double sigma, double TTM, double q, 
		int N, int M, double S_max);
	IFDM(const IFDM& ifdm);
	IFDM& operator=(const IFDM& ifdm);
	~IFDM();

	Vector solve();

private:
	int K;
	double r;
	double sigma;
	double TTM;
	double q;
	
	double delta_tau;
	double delta_x;
	double S_max;

	int N;	// y�� point ���� - 1
	int M;	// x�� point ���� - 1

};

void testIFDM();