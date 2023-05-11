#include "Utils.h"

double simpleDeterminant(const MatrixXd& matrix)
{
	if (matrix.size() == 4) {
		return matrix(0, 0) * matrix(1, 1) -
			matrix(0, 1) * matrix(1, 0);
	}

	double determinant = 0;
	MatrixXd subMatrix;
	for (size_t j = 0; j < matrix.cols(); j++) {
		subMatrix = MatrixXd::Zero(matrix.rows() - 1, matrix.cols() - 1);
		getMinor(matrix, subMatrix, j);
		double postDeterminant = simpleDeterminant(subMatrix);
		determinant += std::powf(-1, 1 + (j + 1)) * matrix(0, j) * postDeterminant;
	}
	return determinant;
}

double threadDetermenant(const MatrixXd& matrix, const uint16_t threadCount)
{
	int htc = boost::thread::hardware_concurrency();
	if (threadCount > htc) {
		std::printf("Превышено кол-во потоков. Выбрано %d из %d", threadCount, htc);
		exit(-1);
	}

	boost::asio::thread_pool pool(threadCount);

	double determinant = 0;

	for (size_t j = 0; j < matrix.cols(); j++)
	{
		MatrixXd subMatrix{ MatrixXd::Zero(matrix.rows() - 1, matrix.cols() - 1) };
		getMinor(matrix, subMatrix, j);

		boost::asio::post(pool, boost::bind(calculateAddition, j, subMatrix, boost::atomic_ref<double>(determinant), matrix(0, j)));
	}
	pool.join();
	return determinant;
}

double mpiDeterminant(const MatrixXd& matrix, int &numprocs, int &rank)
{
	double determinant = 0;
	/*if (rank == 0) {
		int proc = 0;
		for (size_t i = 0; i < matrix.cols(); i++)
		{
			if (proc <= numprocs) {
				MPI_Send(&i, 1, MPI_INTEGER, proc, 1, MPI_COMM_WORLD);
			}
			else {

			}
		}
	}
	else {
		MatrixXd subMatrix{ MatrixXd::Zero(matrix.rows() - 1, matrix.cols() - 1) }
	}*/
	
	return determinant;
}