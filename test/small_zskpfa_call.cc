#include <iostream>
#include <complex>
#include "../src/skpfa.hh"
typedef std::complex<double> dcomplex;

#define matA(i, j) matA[ (i) + (j)*8 ]
#define matC(i, j) matC[ (i) + (j)*8 ]
#define matD(i, j) matD[ (i) + (j)*8 ]
#define mat1(i, j) mat1[ (i) + (j)*8 ]
#define mat2(i, j) mat2[ (i) + (j)*8 ]
#define mat3(i, j) mat3[ (i) + (j)*8 ]

int main(const int argc, const char *argv[])
{
    // these two doesn't conflict.
    using namespace std;

    // prepares input and output container.
    complex<double> matA[8 * 8];
    complex<double> mat1[8 * 4];
    complex<double> mat2[8 * 4];
    complex<double> matC[8 * 8];
    complex<double> matD[8 * 8];
    complex<double> mat3[8 * 4];
    for (signed j = 0; j < 8; ++j)
        for (signed i = 0; i < 8; ++i)
            matA(i, j) = complex<double>(i*2.0 + j, i*0.25 - j*0.5);
    complex<double> alpha = 1.0, beta = 0.0;
    for (unsigned i = 0; i < 8; ++i) {
        for (unsigned j = 0; j < 8; ++j)
            printf("(%10f,%10f) ", real(matA(i, j)), imag(matA(i, j)));
        cout << '\n';
    }

    // call Pfaffian calculation.
    auto Pfa = skpfa<complex<double> >('U', 8,
                               &matA(0, 0), 8, /*inv=*/1,
                               &mat1(0, 0), &mat2(0, 0),
                               &matC(0, 0), &matD(0, 0), &mat3(0, 0), 3);

    printf("Pfa = (%16.8e, %16.8e)\n", real(Pfa), imag(Pfa));
    // row-col formatted output of C.
    for (unsigned i = 0; i < 8; i++) {
        for (unsigned j = 0; j < 8; j++)
            printf("(%10f,%10f) ", real(matA(i, j)), imag(matA(i, j)));
        cout << '\n';
    }
    return 0;
}
