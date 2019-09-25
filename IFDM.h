#pragma once

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

// 일단 european put option으로 작성.
// 나중에 다른 파생상품들도 활용할 수 있도록 확장해보자.

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

	int N;	// y축 point 갯수 - 1
	int M;	// x축 point 갯수 - 1

};

void testIFDM();