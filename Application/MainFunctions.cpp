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

		boost::asio::defer(pool, boost::bind(calculateAddition, j, subMatrix, boost::atomic_ref<double>(determinant), matrix(0, j)));
	}
	pool.join();
	return determinant;
}

double mpiDeterminant(MatrixXd *matrix, int numOfWorkers, int rank, int owner)
{
	double determinant = 0;
	int SEND{ 1 };
	int STOP{ 2 };
	// for owner
	if (!rank) {
		MPI_Status* status = new MPI_Status;;
		//MatrixXd subMatrix;

		int *buff_send_owner = new int(0);
		double *buff_recv_owner = new double(0);

		std::vector<int> clients{0,0,0,1,0};

		for (int i = 0; i < matrix->cols();)
		{
			for (size_t j = 1; j < 5; j++) {
				if (clients[j] != 1) {
					std::cout << "Owner sending:" << i << std::endl;
					MPI_Send(&i, 1, MPI_INTEGER, j, SEND, MPI_COMM_WORLD);
					i++;
					clients[j] = 1;
				}
			}
			
			MPI_Recv(buff_recv_owner, 1, MPI_INTEGER, MPI_ANY_SOURCE, SEND, MPI_COMM_WORLD, status);
			std::printf("Owner recevied %d data from %d", *buff_recv_owner, status->MPI_SOURCE);
			clients[status->MPI_SOURCE] = 0;
		}
		//	//getMinor(*matrix, subMatrix, i);
		//	determinant += std::powf(-1, 1 + (i + 1)) * (*matrix)(0, i) * simpleDeterminant(subMatrix);
		//	i++;
		//	std::cout << "owner  listen" << rank;
		//	MPI_Recv(buff_recv_owner, 1, MPI_DOUBLE, MPI_ANY_SOURCE, SEND, MPI_COMM_WORLD, status);
		//	clients[status->MPI_SOURCE] = 0;
		//	determinant += *buff_recv_owner;
		//}
		//delete buff_recv_owner, buff_send_owner, status;
		//for (size_t j = 1; j < clients.size(); j++) {
		//		std::cout << "owner  send" << rank;
		//		*buff_send_owner = j;
		//		MPI_Send(buff_send_owner, 1, MPI_INTEGER, j, STOP, MPI_COMM_WORLD);
		//}
	}

	
	// for workers
	if (rank) {
		//MPI_Status* status = new MPI_Status;
		//MatrixXd subMatrix;
		//double* buff_send_worker = new double(0);
		//int* buff_recv_worker = new int(0);
		//while (true) {
		//	std::cout << "client  listen" << rank;

		//	MPI_Recv(buff_recv_worker, 1, MPI_INTEGER, owner, SEND, MPI_COMM_WORLD, status);
		//	if (status->MPI_TAG == STOP) { break; }

		//	//getMinor(*matrix, subMatrix, *buff_recv_worker);
		//	*buff_send_worker = std::powf(-1, 1 + (*buff_recv_worker + 1)) * (*matrix)(0, *buff_recv_worker) * simpleDeterminant(subMatrix);
		//	std::cout << "client  send" << rank;
		//	MPI_Send(buff_send_worker, 1, MPI_DOUBLE, owner, SEND, MPI_COMM_WORLD);
		//}
		//delete status, buff_send_worker, buff_recv_worker;
	/*	while (true) {
			MPI_Status status;
			int t = 0;
			
			MPI_Recv(&t, 1, MPI_INTEGER, owner, SEND, MPI_COMM_WORLD, &status);
			MatrixXd subMatrix{ Matrix::Ones((matrix->rows() - 1, matrix->cols() - 1)) };
			getMinor(*matrix, subMatrix, t);
			std::cout << subMatrix << std::endl;
			std::printf("Worker %d receved data", t);
			MPI_Send(&t, 1, MPI_INTEGER, owner, SEND, MPI_COMM_WORLD);
			std::printf("Worker %d sended data", t);
		}*/
	}
	return determinant;
}