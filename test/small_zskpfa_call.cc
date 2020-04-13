#include <iostream>
#include <complex>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "../src/skpfa.hh"

int main(const int argc, const char *argv[])
{
    // these two doesn't conflict.
    using namespace std;
    using namespace boost::numeric::ublas;

    // prepares input and output container.
    matrix<complex<double>, column_major> matA(8, 8);
    matrix<complex<double>, column_major> mat1(8, 4);
    matrix<complex<double>, column_major> mat2(8, 4);
    matrix<complex<double>, column_major> matC(8, 8);
    matrix<complex<double>, column_major> matD(8, 8);
    matrix<complex<double>, column_major> mat3(8, 4);
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

    printf("Pfa = %16.8e\n", real(Pfa));
    // row-col formatted output of C.
    for (unsigned i = 0; i < 8; i++) {
        for (unsigned j = 0; j < 8; j++)
            printf("(%10f,%10f) ", real(matA(i, j)), imag(matA(i, j)));
        cout << '\n';
    }
    return 0;
}
