#pragma once
#include <Eigen/Dense>
#include <iostream>
#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>
#include <boost/thread.hpp>
#include <boost/atomic/atomic_ref.hpp>
#include <boost/timer/timer.hpp>
#include <mpi.h>
#include <string>

using namespace Eigen;



double simpleDeterminant(const MatrixXd& matrix);
double threadDetermenant(const MatrixXd& matrix, const uint16_t threadCount);
double mpiDeterminant(const MatrixXd& matrix, int& numrocs, int& rank);
void getMinor(const MatrixXd& matrix, MatrixXd& subMatrix, const uint16_t j);
void calculateAddition(uint16_t j, MatrixXd subMatrix, boost::atomic_ref<double> determinant, double point);

