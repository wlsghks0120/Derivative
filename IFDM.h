#pragma once

#include "FDMutils.h"
#include <vector>

// �ϴ� european put option���� �ۼ�.
// ���߿� �ٸ� �Ļ���ǰ�鵵 Ȱ���� �� �ֵ��� Ȯ���غ���.

class IFDM {
public:
	IFDM();
	IFDM(int K, double r, double sigma, double TTM, double q, 
		int N, int M, double S_max);
	IFDM(const IFDM& ifdm);
	IFDM& operator=(const IFDM& ifdm);
	~IFDM();

	std::unique_ptr<Vector> solve();

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