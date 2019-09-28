#pragma once

#include <vector>
#include "FDMutils.h"

/*
���� ������ġ�� �ȸ��������. std::vector �޴� �����ڿ���
A�� d�� ����� ��ġ�ϴ��� Ȯ���� ���ϰ� �ְ�
����� Tridiagonal ����, Square����, Strictly diagonally dominant ���� Ȯ�� �ؾ���.
*/

class ThomasAlgorithm {
public:
	ThomasAlgorithm();
	ThomasAlgorithm(const std::unique_ptr<Matrix>& A, const std::unique_ptr<Vector>& d);
	~ThomasAlgorithm();

	int eqtSize();
	std::unique_ptr<Vector> solve();

private:
	std::unique_ptr<Matrix> A;
	std::unique_ptr<Vector> d;

	// ���� �̱��� �� ������ġ �Լ���
	bool isSquare();
	bool isTridiagonal();
	bool isStrictlyDiagonallyDominant();
};

void testThomasAlgorithm();
