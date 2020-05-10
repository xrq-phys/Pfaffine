#include <cstdio>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include "../src/skr2k.tcc"
#include "../src/blalink.hh"

#define M_ 64
#define N_ 60
#define K_ 8

int main(const int argc, const char *argv[])
{
    using namespace boost::numeric::ublas;
    using namespace boost::random;
    // Plain C IO.
    auto *fid_a = fopen("AoT.dat", "w");
    auto *fid_b = fopen("BoT.dat", "w");
    auto *fid_c1 = fopen("Co1.dat", "w");
    auto *fid_c2 = fopen("Co2.dat", "w");

    // Initialize random.
    mt19937 rng;
    normal_distribution<double> randn(0.0, 1.0);

    matrix<double, column_major> matA(M_, K_);
    matrix<double, column_major> matB(N_, K_);
    matrix<double, column_major> matC(M_, N_);
    double *ptrA = &(matA(0, 0));
    double *ptrB = &(matB(0, 0));
    double *ptrC = &(matC(0, 0));
    double blasp[700];
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

#ifdef _SVE
    printf("Skipping plain.\n");
#else
    printf("Executing plain.\n");
    for (unsigned j = 0; j < N_; ++j)
        for (unsigned i = 0; i < M_; ++i)
            matC(i, j) = 0.0;
    gemm('N', 'T', M_, N_, K_, 1.0,
         ptrA, M_,
         ptrB, N_, 1.0,
         ptrC, M_);
    for (unsigned i = 0; i < M_; ++i) {
        for (unsigned j = 0; j < N_; ++j)
            fprintf(fid_c1, "%16.8e ", matC(i, j));
        fprintf(fid_c1, "\n");
    }
#endif

    // Clear C and apply skr2k.
    for (unsigned j = 0; j < N_; ++j)
        for (unsigned i = 0; i < M_; ++i)
            matC(i, j) = 0.0;
    unsigned mr, nr;
    set_blk_size<double>(&mr, &nr);
    printf("Executing custom.\n");
    ogemm<double>(M_, N_, K_, 1.0, 
                  ptrA, M_, 
                  ptrB, N_, 1.0, 
                  ptrC, M_, mr, nr, blasp);
    for (unsigned i = 0; i < M_; ++i) {
        for (unsigned j = 0; j < N_; ++j)
            fprintf(fid_c2, "%16.8e ", matC(i, j));
        fprintf(fid_c2, "\n");
    }

    fclose(fid_a);
    fclose(fid_b);
    fclose(fid_c1);
    fclose(fid_c2);
    return 0;
}
