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
    std::vector<double> simpleAlgTimeForThread(0);
    std::vector<double> simpleAlgVarianceInForThread(0);

    
    std::vector<double> threadAlgTime(0);
    std::vector<double> threadAlgTimeInForThread(0);
    std::vector<double> threadAlgVarianceInForThread(0);


    for (size_t thread_s = 1; thread_s <= thread_max; thread_s++)
    {
        std::printf("Thread number:[%zd]:\n\n", thread_s);
        threadAlgTime.clear();
        for (size_t sampl = 0; sampl < sample_max; sampl++)
        {
            MatrixXd matrix{ 10 * MatrixXd::Random(n,n) };
            timer.start();
            double resultSimple{ simpleDeterminant(matrix) };
            timer.stop();

            simpleAlgTime.push_back(atof(timer.format(6, "%w").c_str()));
            std::printf("\n\tSimpleTime: %.3f s\n", simpleAlgTime.at(sampl));

            timer.start();
            double resultThread{ threadDetermenant(matrix, thread_s) };
            timer.stop();

            threadAlgTime.push_back(atof(timer.format(6, "%w").c_str()));
            std::printf("\n\tthreadTime: %.3f s\n", threadAlgTime.at(sampl));
        }

        
        accumulator_set<double, features<tag::mean, tag::variance>> accThreadMean;
        accumulator_set<double, features<tag::mean, tag::variance>> accSimpleMean;
        accumulator_set<double, features<tag::mean>> accMPIMean;

        for (auto var : threadAlgTime) {
            accThreadMean(var);
        }
        
        for (auto var : simpleAlgTime) {
            accSimpleMean(var);
        }

        threadAlgTimeInForThread.push_back(boost::accumulators::mean(accThreadMean));
        threadAlgVarianceInForThread.push_back(boost::accumulators::variance(accThreadMean));
        
        simpleAlgTimeForThread.push_back(boost::accumulators::mean(accSimpleMean));
        simpleAlgVarianceInForThread.push_back(boost::accumulators::variance(accSimpleMean));


        std::printf("\nMean time thread:%.3f\n\n", boost::accumulators::mean(accThreadMean));
        std::printf("\nMean time simple:%.3f\n\n", boost::accumulators::mean(accSimpleMean));
        //for () {}

    }

    //accumulator_set<double, features<tag::mean, tag::variance>> accMPI;


    //double mpiMean{ 0 };

    //double mpiVariance{ 0 };

    ///for () {}
    
    students_t dist(thread_max - 1);
    double T = boost::math::quantile(complement(dist, 0.05 / 2));
    double wSimple = T * simpleAlgTimeForThread.at(thread_max - 1) / std::sqrt(double(simpleAlgTimeForThread.size()));
    double wThread = T * threadAlgVarianceInForThread.at(thread_max - 1) / std::sqrt(double(threadAlgTimeInForThread.size()));

    //double wMPI{ 0 };

    double lowerLimitSimple = simpleAlgTimeForThread.at(thread_max - 1) - wSimple;
    double upperLimitSimple = simpleAlgTimeForThread.at(thread_max - 1) + wSimple;
    
    double lowerLimitThread = threadAlgVarianceInForThread.at(thread_max - 1) - wThread;
    double upperLimitThread = threadAlgVarianceInForThread.at(thread_max - 1) + wThread;
    
    /*double lowerLimitThread = threadMean - wThread;
    double upperLimitThread = threadMean + wThread;*/


    std::printf("\t\t\t\tAll statistic:\n\n\nthread_number\tthread_mean\tsimple_mean\tmpi_mean\n\n");
    for (int i = 1; i <= thread_max; i++)
    {
        std::printf("\t[%d]\t\t%.3e\t\t%.3f\t%.3e\n", i, threadAlgTimeInForThread.at(i -1 ), simpleAlgTimeForThread.at(i - 1), 0.0);
    }
    std::printf("Statistic for %d samples, p:%.2f%%:\n", sample_max, 95.);
    std::printf("\tThread: [mean: %.3e, variance: %.3e,confidence interval: (%.3e;%.3e)]\n",threadAlgTimeInForThread.at(thread_max -1),
        threadAlgVarianceInForThread.at(thread_max - 1), lowerLimitThread, upperLimitThread);
    std::printf("\tSimple: [mean: %.3e, variance: %.3e,confidence interval: (%.3e;%.3e)]\n", simpleAlgTimeForThread.at(thread_max - 1),
        simpleAlgVarianceInForThread.at(thread_max - 1),
        lowerLimitSimple, upperLimitSimple);

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

