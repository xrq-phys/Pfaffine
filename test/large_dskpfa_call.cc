#include <iostream>
#include <random>
#include <chrono>
#include "../src/skpfa.hh"
#ifdef _Intel_Advisor
#include <ittnotify.h>
#endif
#ifdef _Fujitsu_Profiler
#include <fj_tool/fapp.h>
#endif

static const unsigned N = 2240;
static const unsigned NPANEL = 16;

#define matA(i, j) matA[ (i) + (j)*N ]
#define mat1(i, j) mat1[ (i) + (j)*N ]
#define matC(i, j) matC[ (i) + (j)*N ]

int main(const int argc, const char *argv[])
{
    // these two doesn't conflict.
    using namespace std;

    // use random.
    mt19937_64 rng(1234);
    normal_distribution<double> dist(0, 0.037);

#ifdef _Intel_Advisor
    // initialization phase.
    // stop data collection.
    __itt_pause();
#endif

    // prepares input and output container.
    double matA[N * N];
    double mat1[N * NPANEL];
    double matC[N * N];
    signed iPov[N + 1];
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

#ifdef _Intel_Advisor
    // start Advisor's collection.
    __itt_resume();
#endif
#ifdef _Fujitsu_Profiler
    fapp_start("skpfa", 1, 0);
#endif

    // call Pfaffian calculation.
    auto start = std::chrono::high_resolution_clock::now();
    double Pfa;
    signed info = skpfa<double>(BLIS_UPPER, N,
                                &matA(0, 0), N,
                                &matC(0, 0), N,
                                iPov, false, &Pfa,
                                &mat1(0, 0), N * NPANEL);
    auto elapsed = std::chrono::high_resolution_clock::now() - start;
    long long microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();

#ifdef _Intel_Advisor
    // stop collection.
    __itt_pause();
#endif
#ifdef _Fujitsu_Profiler
    fapp_stop("skpfa", 1, 0);
#endif

    printf("Pfa = %16.8e\n", Pfa);
    printf("Time consumed: %20lf ms\n", double(microseconds)/1000);
    // row-col formatted output of C.
    /*
    for (unsigned i = 0; i < 8; i++) {
        for (unsigned j = 0; j < 8; j++)
            printf("%12f ", matA(i, j));
        cout << '\n';
    }*/
    return 0;
}
