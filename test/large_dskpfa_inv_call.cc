#include <iostream>
#include <random>
#include "../src/skpfa.hh"

static const unsigned N = 128;
static const unsigned NPANEL = 8;

#define matA(i, j) matA[ (i) + (j)*N ]
#define mat1(i, j) mat1[ (i) + (j)*N ]
#define mat2(i, j) mat2[ (i) + (j)*N ]
#define matC(i, j) matC[ (i) + (j)*N ]

int main(const int argc, const char *argv[])
{
    // these two doesn't conflict.
    using namespace std;

    // use random.
    mt19937_64 rng(1234);
    normal_distribution<double> dist(0, 0.05);

    // prepares input and output container.
    double matA[N * N];
    double mat1[N * NPANEL];
    double mat2[N * NPANEL];
    double matC[N * N];
    double matD[N * N];
    double mat3[N * (NPANEL+1)];
    for (unsigned j = 0; j < N; ++j)
        for (unsigned i = 0; i < N; ++i)
            matA(i, j) = dist(rng);
    double alpha = 1.0, beta = 0.0;
    // Write input to file.
    // As N*N is too large.
    char fnm[20];
    sprintf(fnm, "IN%d.dat", N);
    auto *fin = fopen(fnm, "w");
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned j = 0; j < N; ++j)
            fprintf(fin, "%12f ", matA(i, j));
        fprintf(fin, "\n");
    }
    fclose(fin);

    // call Pfaffian calculation.
    double Pfa = skpfa<double>('U', N,
                               &matA(0, 0), N, /*inv=*/1,
                               &mat1(0, 0), &mat2(0, 0),
                               //(nullptr, nullptr, nullptr, NPANEL);
                               matC, matD, mat3, NPANEL);

    printf("Pfa = %16.8e\n", Pfa);
    // row-col formatted output of C.
    sprintf(fnm, "OUT%d.dat", N);
    auto *fout = fopen(fnm, "w");
    for (unsigned i = 0; i < N; i++) {
        for (unsigned j = 0; j < N; j++)
            fprintf(fout, "%12f ", matA(i, j));
        fprintf(fout, "\n");
    }
    fclose(fout);

    return 0;
}
