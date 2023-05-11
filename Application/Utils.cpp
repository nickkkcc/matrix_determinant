#include "Utils.h"

void calculateAddition(uint16_t j, MatrixXd subMatrix, boost::atomic_ref<double> determinant, double point) {
	determinant += std::powf(-1, 1 + (j + 1)) * point* simpleDeterminant(subMatrix);
}

void getMinor(const MatrixXd& matrix, MatrixXd& subMatrix, uint16_t j)
{
	if (j == 0) {
		subMatrix = matrix.block(1, 1, matrix.rows() - 1, matrix.cols() - 1);
	}
	else if (j == matrix.cols() - 1) {

		subMatrix = matrix.block(1, 0, matrix.rows() - 1, matrix.cols() - 1);
	}
	else {
		subMatrix << matrix.block(1, 0, matrix.rows() - 1, j),
			matrix.block(1, j + 1, matrix.rows() - 1, matrix.cols() - j - 1);
	}
}



