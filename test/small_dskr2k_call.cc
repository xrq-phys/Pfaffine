#include <cstdio>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include "../src/skr2k.hh"
#include "../src/blalink.hh"

#define M_ 64
#define N_ 64
#define K_ 8

int main(const int argc, const char *argv[])
{
    using namespace boost::numeric::ublas;
    using namespace boost::random;
    // Plain C IO.
    auto *fid_a = fopen("AT.dat", "w");
    auto *fid_b = fopen("BT.dat", "w");
    auto *fid_c1 = fopen("C1.dat", "w");
    auto *fid_c2 = fopen("C2.dat", "w");

    // Initialize random.
    mt19937 rng;
    normal_distribution<double> randn(0.0, 1.0);

    matrix<double, column_major> matA(M_, K_);
    matrix<double, column_major> matB(N_, K_);
    matrix<double, column_major> matC(M_, N_);
    double *ptrA = &(matA(0, 0));
    double *ptrB = &(matB(0, 0));
    double *ptrC = &(matC(0, 0));
    for (unsigned l = 0; l < K_; ++l) {
        for (unsigned i = 0; i < M_; ++i) {
            matA(i, l) = randn(rng);
            fprintf(fid_a, "%16.8e ", matA(i, l));
        } fprintf(fid_a, "\n");
    }
    for (unsigned l = 0; l < K_; ++l) {
        for (unsigned j = 0; j < N_; ++j) {
            matB(j, l) = randn(rng);
            fprintf(fid_b, "%16.8e ", matB(j, l));
        } fprintf(fid_b, "\n");
    }

    printf("Executing plain.\n");
    for (unsigned j = 0; j < N_; ++j)
        for (unsigned i = 0; i < M_; ++i)
            matC(i, j) = 0.0;
    gemm('N', 'T', M_, N_, K_, 1.0,
         ptrA, M_,
         ptrB, N_, 1.0,
         ptrC, M_);
    gemm('N', 'T', N_, M_, K_, -1.0,
         ptrB, N_,
         ptrA, M_, 1.0,
         ptrC, M_);
    for (unsigned i = 0; i < N_; ++i) {
        for (unsigned j = 0; j < M_; ++j)
            fprintf(fid_c1, "%16.8e ", matC(i, j));
        fprintf(fid_c1, "\n");
    }

    // Clear C and apply skr2k.
    for (unsigned j = 0; j < N_; ++j)
        for (unsigned i = 0; i < M_; ++i)
            matC(i, j) = 0.0;
    printf("Executing accelerated.\n");
    skr2k<double>('U', 'N', N_, K_, 1.0, 
                  ptrA, M_, 
                  ptrB, N_, 1.0, 
                  ptrC, M_);
    for (unsigned i = 0; i < N_; ++i) {
        for (unsigned j = 0; j < M_; ++j)
            fprintf(fid_c2, "%16.8e ", matC(i, j));
        fprintf(fid_c2, "\n");
    }

    fclose(fid_a);
    fclose(fid_b);
    fclose(fid_c1);
    fclose(fid_c2);
    return 0;
}

