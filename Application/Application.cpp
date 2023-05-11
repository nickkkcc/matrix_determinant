#include <iostream>
#include "Utils.h"
#include <vector>
#include <boost/math/distributions/students_t.hpp>
#include "boost/accumulators/accumulators.hpp"
#include "boost/accumulators/statistics.hpp"


using namespace boost::accumulators;
using namespace boost::math;

void validate(int argc, char** argv) {
    if (strcmp("-matrix_size", argv[1]) <  0 || strcmp("-max_thread_count", argv[3]) < 0 || atoi(argv[2]) == 0 || atoi(argv[4]) == 0 ||
        strcmp("-samp_size", argv[5]) < 0 || atoi(argv[6]) == 0) {
        std::cout << "Wrong arguments!\n Write mpiexec -n [numberProc] program_name.exe -matrix_size [size]{matrix: size X size} -max_thread_count [n] -samp_size [size]";
        exit(-1);
    }
}


int main(int argc, char** argv)
{
    //int rank, numprocs;
    //MPI_Init(&argc, &argv);
    //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    // -------------------------------------------------

    validate(argc, argv);
    int n{ atoi(argv[2]) };
    int thread_max{ atoi(argv[4]) };
    int sample_max{ atoi(argv[6]) };

    boost::timer::cpu_timer timer;

    std::vector<double> simpleAlgTime(0);

    
    std::vector<double> threadAlgTime(0);
    std::vector<double> threadAlgTimeInForThread(0);


    for (size_t thread_s = 1; thread_s <= thread_max; thread_s++)
    {
        std::printf("Thread number:[%zd]:\n\n", thread_s);
        threadAlgTime.clear();
        for (size_t sampl = 0; sampl < sample_max; sampl++)
        {
            MatrixXd matrix{ 10 * MatrixXd::Random(n,n) };
            if (thread_s == thread_max) {
                
                timer.start();
                double resultSimple{ simpleDeterminant(matrix) };
                timer.stop();

                simpleAlgTime.push_back(atof(timer.format(6, "%w").c_str()));
                std::printf("\n\tSimpleTime: %.3e s\n", simpleAlgTime.at(sampl));
            }

            timer.start();
            double resultThread{ threadDetermenant(matrix, thread_s) };
            timer.stop();

            threadAlgTime.push_back(atof(timer.format(6, "%w").c_str()));
            std::printf("\n\tthreadTime: %.3e s\n", threadAlgTime.at(sampl));
        }

        
        accumulator_set<double, features<tag::mean>> accThreadMean;
        accumulator_set<double, features<tag::mean>> accMPIMean;

        for (auto var : threadAlgTime) {
            accThreadMean(var);
        }
        threadAlgTimeInForThread.push_back(boost::accumulators::mean(accThreadMean));
        std::printf("\nMean time thread:%.3e\n\n", boost::accumulators::mean(accThreadMean));
        //for () {}

    }

    accumulator_set<double, features<tag::mean, tag::variance>> accSimple;
    accumulator_set<double, features<tag::mean, tag::variance>> accThread;
    //accumulator_set<double, features<tag::mean, tag::variance>> accMPI;

    for (auto var : threadAlgTimeInForThread) {
        accThread(var);
    }
    for (auto var : simpleAlgTime) {
        accSimple(var);
    }

    double simpleMean{ boost::accumulators::mean(accSimple) };
    double threadMean {boost::accumulators::mean(accThread)};
    //double mpiMean{ 0 };

    double simpleVariance{ boost::accumulators::variance(accSimple) };
    double threadVariance{ boost::accumulators::variance(accThread) };
    //double mpiVariance{ 0 };

    ///for () {}
    
    students_t dist(thread_max - 1);
    double T = boost::math::quantile(complement(dist, 0.05 / 2));
    double wSimple = T * simpleVariance / std::sqrt(double(simpleAlgTime.size()));
    double wThread = T * threadVariance / std::sqrt(double(threadAlgTimeInForThread.size()));

    //double wMPI{ 0 };

    double lowerLimitSimple = simpleMean - wSimple;
    double upperLimitSimple = simpleMean + wSimple;
    
    double lowerLimitThread = threadMean - wThread;
    double upperLimitThread = threadMean + wThread;
    
    /*double lowerLimitThread = threadMean - wThread;
    double upperLimitThread = threadMean + wThread;*/


    std::printf("\t\t\t\tAll statistic:\n\n\nthread_number\tsimple_mean\tthread_mean\tmpi_mean\n\n");
    for (int i = 1; i <= thread_max; i++)
    {
        std::printf("\t[%d]\t\t%.3f\t\t%.3f\t%.3f\n", i, threadAlgTime.at(i), simpleMean, 0.0);
    }
    std::printf("Statistic for %d samples, p:%.2f%%:\n", sample_max, 95.);
    std::printf("\tThread: [mean: %.3e, variance: %.3e,confidence interval: (%.3e;%.3e)]\n",threadMean, threadVariance, lowerLimitThread, upperLimitThread);
    std::printf("\tSimple: [mean: %.3e, variance: %.3e,confidence interval: (%.3e;%.3e)]\n",simpleMean, simpleVariance, lowerLimitSimple, upperLimitSimple);

    // -------------------------------------------------

    //MPI_Status status;

    //// Execute thread and non-thread alghorims for 1 process.
    //if (rank == 0) {

     // -------------------------------------------------

      

        // -------------------------------------------------
        //int flag = 1;
        //for (size_t i = 1; i < numprocs; i++)
        //{

        //    MPI_Send(&flag,1,MPI_INTEGER,i,0,MPI_COMM_WORLD);
        //}
        //
    //}
    //// For other
    //else {
    //    int flag;
    //    MPI_Recv(&flag, 1, MPI_INTEGER, 0, 0,MPI_COMM_WORLD, &status);
    //}
    //if (rank == 0) {
    //    timer.start();
    //}
    //double resultMPI{ mpiDeterminant(matrix, numprocs, rank) };
    //if (rank == 0) {
    //    timer.stop();
    //    std::cout << timer.format(3, "mpiTime: %w sec\n");
    //}
    //MPI_Finalize();
    return 0;
}

