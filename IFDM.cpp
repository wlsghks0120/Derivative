#include "IFDM.h"
#include "ThomasAlgorithm.h"
#include "TridiagonalUtils.h"

#include <time.h>
#include <cmath>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

using std::cout;
using std::endl;

void testIFDM()
{
	int K = 10;
	double r = 0.05;
	double sigma = 0.2;
	double TTM = 1. / 2.;
	double q = 0;
	int N = 1000;
	int M = 10000;
	double S_max = 20;
	
	clock_t start, end;
	start = clock();
	IFDM ifdm(K, r, sigma, TTM, q, N, M, S_max);
	IFDM::Vector price = ifdm.solve();
	end = clock();
	cout << "option price is : " << price(N / 2) << endl;
	cout << "time = " << double(end-start) << "ms" << endl;
}

IFDM::IFDM()
	: K(0), r(0), sigma(0), TTM(0), q(0), delta_tau(0), delta_x(0), S_max(0), N(0), M(0)
{
}

IFDM::IFDM(int K, double r, double sigma, double TTM, double q,
	int N, int M, double S_max)
	: K(K), r(r), sigma(sigma), TTM(TTM), q(q), S_max(S_max), N(N), M(M), 
	delta_x(S_max/N), delta_tau(TTM/M)
{
}

IFDM::IFDM(const IFDM& ifdm)
	: K(ifdm.K), r(ifdm.r), sigma(ifdm.sigma), TTM(ifdm.TTM), q(ifdm.q), 
	delta_tau(ifdm.delta_tau), delta_x(ifdm.delta_x), S_max(ifdm.S_max), N(ifdm.N), M(ifdm.M)
{
}

IFDM& IFDM::operator=(const IFDM& ifdm)
{
	if (this != &ifdm) {
		K = ifdm.K;
		r = ifdm.r;
		sigma = ifdm.sigma;
		TTM = ifdm.TTM;
		q = ifdm.q;
		delta_tau = ifdm.delta_tau;
		delta_x = ifdm.delta_x;
		S_max = ifdm.S_max;
		N = ifdm.N;
		M = ifdm.M;
	}
	return *this;
}

IFDM::~IFDM()
{
}

IFDM::Vector IFDM::solve()
{
	cout << "IFDM::solve() is called" << endl;
	// put option payoff 사용.
	// x축, y축 point 갯수 계산
	Matrix v(N+1, M+1);
	for (int n = 0; n <= N; ++n) {
		v(n, 0) = std::max(K - (N - n) * delta_x, 0.0);	// 격자 왼쪽 벽 위부터 채우기
	}
	for (int m = 0; m <= M; ++m) {
		v(0, m) = 0;		// 격자 위쪽 벽 왼쪽부터 채우기
		v(N, m) = K * std::exp(-r * delta_tau * m) - 0;		// 격자 아래쪽 벽 왼쪽부터 채우기
	}

	// Ax = d를 풀기 위해 A,d를 셋팅. a,b,c,는 A의 대각 벡터들
	Vector d(N+1);
	for (int i = 0; i <= N; ++i)
		d(i) = v(N-i, 0);
	Vector a(N+1), b(N), c(N);	
	
	// a의 양끝은 매 step마다 따로 계산
	a(0) = 0;
	a(N) = 1;	// 풋옵션이라 이렇게. 맨 윗줄이 다 0이어서.
	c(0) = 0;
	for (int i = 1; i <= N - 1; ++i) {
		a(i) = 1. + delta_tau * (r + sigma * sigma * i * i);
		c(i) = delta_tau / 2. * (-(r - q) * i - sigma * sigma * i * i);
	}
	b(N - 1) = 0;
	for (int i = 0; i <= N - 2; ++i) {
		b(i) = delta_tau / 2 * ((r - q) * (i + 1) - sigma * sigma * (i + 1) * (i + 1));
	}
	Matrix A = createTridiagonalMatrix(a, b, c);

	for (int m = 1; m <= M; ++m) {
		A(0, 0) = v(N, m - 1) / v(N, m);
		if (m % 200 == 0) cout << m << "th step is completed" << endl;
		ThomasAlgorithm ta(A, d);
		d = ta.solve();
		for (int i = 0; i <= N; ++i)
			v(N - i,m) = d(i);
	}

	return d;
}

