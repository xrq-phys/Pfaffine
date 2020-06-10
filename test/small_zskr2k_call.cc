#include <cstdio>
#include <complex>
#include <random>
#include "../src/skr2k.hh"
#include "../src/blalink.hh"
typedef std::complex<double> dcomplex;

#define M_ 64
#define N_ 64
#define K_ 8

#define matA(i, j) matA[ (i) + (j)*M_ ]
#define matB(i, j) matB[ (i) + (j)*N_ ]
#define matC(i, j) matC[ (i) + (j)*M_ ]

int main(const int argc, const char *argv[])
{
    using namespace std;
    // Plain C IO.
    auto *fid_a = fopen("AT.dat", "w");
    auto *fid_b = fopen("BT.dat", "w");
    auto *fid_c1 = fopen("C1.dat", "w");
    auto *fid_c2 = fopen("C2.dat", "w");

    // Initialize random.
    mt19937 rng;
    normal_distribution<double> randn(0.0, 1.0);

    dcomplex matA[M_ * K_];
    dcomplex matB[N_ * K_];
    dcomplex matC[M_ * N_];
    dcomplex *ptrA = &(matA(0, 0));
    dcomplex *ptrB = &(matB(0, 0));
    dcomplex *ptrC = &(matC(0, 0));
    dcomplex blasp[600];
    for (unsigned l = 0; l < K_; ++l) {
        for (unsigned i = 0; i < M_; ++i) {
            matA(i, l) = dcomplex(randn(rng), randn(rng));
            fprintf(fid_a, "%16.8e %16.8e ", std::real(matA(i, l)), std::imag(matA(i, l)));
        } fprintf(fid_a, "\n");
    }
    for (unsigned l = 0; l < K_; ++l) {
        for (unsigned j = 0; j < N_; ++j) {
            matB(j, l) = dcomplex(randn(rng), randn(rng));
            fprintf(fid_b, "%16.8e %16.8e ", std::real(matB(j, l)), std::imag(matB(j, l)));
        } fprintf(fid_b, "\n");
    }

    printf("Executing plain.\n");
    for (unsigned j = 0; j < N_; ++j)
        for (unsigned i = 0; i < M_; ++i)
            matC(i, j) = 0.0;
    gemm('N', 'T', M_, N_, K_, dcomplex(1.0),
         ptrA, M_,
         ptrB, N_, dcomplex(1.0),
         ptrC, M_);
    gemm('N', 'T', N_, M_, K_, dcomplex(-1.0),
         ptrB, N_,
         ptrA, M_, dcomplex(1.0),
         ptrC, M_);
    for (unsigned i = 0; i < M_; ++i) {
        for (unsigned j = 0; j < N_; ++j)
            fprintf(fid_c1, "%16.8e %16.8e ", std::real(matC(i, j)), std::imag(matC(i, j)));
        fprintf(fid_c1, "\n");
    }

    // Clear C and apply skr2k.
    for (unsigned j = 0; j < N_; ++j)
        for (unsigned i = 0; i < M_; ++i)
            matC(i, j) = 0.0;
    printf("Executing accelerated.\n");
    skr2k<dcomplex>('U', 'N', N_, K_, dcomplex(1.0),
                    ptrA, M_,
                    ptrB, N_, dcomplex(1.0),
                    ptrC, M_, blasp);
    for (unsigned i = 0; i < M_; ++i) {
        for (unsigned j = 0; j < N_; ++j)
            fprintf(fid_c2, "%16.8e %16.8e ", std::real(matC(i, j)), std::imag(matC(i, j)));
        fprintf(fid_c2, "\n");
    }

    fclose(fid_a);
    fclose(fid_b);
    fclose(fid_c1);
    fclose(fid_c2);
    return 0;
}

