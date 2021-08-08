/*
 * This test has additional dependency:
 * - MPI
 *
 * This test is not compiled by default.
 */
#include <iostream>
#include <random>
#include <chrono>
#include <mpi.h>
#include "../src/skpfa.hh"

const unsigned N_Max = 1200;
const unsigned N_Lst[14] = {
     32,  64,  128,  192, 256,
    320, 384,  448,  512, 640,
    768, 896, 1000, 1200
};
static const unsigned NPANEL = 16;

#define matA(i, j) matA[ (i) + (j)*N ]
#define mat1(i, j) mat1[ (i) + (j)*N ]
#define matC(i, j) matC[ (i) + (j)*N ]

int main(int argc, char *argv[])
{
    using namespace std;
    // MPI information
    int mpi_size, mpi_rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    // use random.
    mt19937_64 rng(1234);
    normal_distribution<double> dist(0, 0.05);

    // prepares input and output container.
    double matA[N_Max * N_Max];
    double mat1[N_Max * NPANEL];
    double matC[N_Max * N_Max];
    signed iPiv[N_Max + 1];

    for (unsigned idx = 0; idx < 14; ++idx) {
        // NOTE: this sets size as well as ldA above.
        int N = N_Lst[idx];

        for (unsigned j = 0; j < N; ++j)
            for (unsigned i = 0; i < N; ++i)
                matA(i, j) = dist(rng);
        double alpha = 1.0, beta = 0.0;
        // Write input to file.
        // As N*N is too large.
        if (N < 600 && mpi_rank == 0) {
            char fnm[20];
            sprintf(fnm, "IN%d.dat", N);
            auto *fin = fopen(fnm, "w");
            for (unsigned i = 0; i < N; ++i) {
                for (unsigned j = 0; j < N; ++j)
                    fprintf(fin, "%12f ", matA(i, j));
                fprintf(fin, "\n");
            }
        }

        auto start = std::chrono::high_resolution_clock::now();
        // call Pfaffian calculation.
        double Pfa;
        signed info = skpfa<double>(BLIS_UPPER, N,
                                    &matA(0, 0), N,
                                    &matC(0, 0), N,
                                    iPiv, true, &Pfa,
                                    &mat1(0, 0), N * NPANEL);
        // collect elapsed time information.
        auto elapsed = std::chrono::high_resolution_clock::now() - start;
        long long microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();

        printf("Rank %3d @ N = %8d: Pfa = %16.8e\n", mpi_rank, N, Pfa);
        printf("Rank %3d @ N = %8d: Time consumed: %20lf ms\n", mpi_rank, N, double(microseconds)/1000);
        fflush(stdout);
    }

    MPI_Finalize();
    return 0;
}
