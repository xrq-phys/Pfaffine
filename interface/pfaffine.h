/**
 * \file pfaffine.h
 * C99/Fortran Interface for Pfaffine calculation.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#ifdef __cplusplus
#define EXTERNC extern "C"
#include <complex>
typedef std::complex<float>  scomplex;
typedef std::complex<double> dcomplex;
#else
#define EXTERNC
#include <complex.h>
typedef float  complex scomplex;
typedef double complex dcomplex;
#endif

// Simplified interface.
EXTERNC void sskpfa(float    *Pfa, char uplo, unsigned n, float    *A, unsigned ldA, unsigned inv);
EXTERNC void dskpfa(double   *Pfa, char uplo, unsigned n, double   *A, unsigned ldA, unsigned inv);
EXTERNC void cskpfa(scomplex *Pfa, char uplo, unsigned n, scomplex *A, unsigned ldA, unsigned inv);
EXTERNC void zskpfa(dcomplex *Pfa, char uplo, unsigned n, dcomplex *A, unsigned ldA, unsigned inv);

// Full interface, callable from both C and Fortran.
EXTERNC void cal_sskpfa_(float    *Pfa, char *uplo, unsigned *n, float    *A, unsigned *ldA, unsigned *inv,
        float    *Sp1, float    *Sp2, float    *Sp3, float    *Sp4, float    *Sp5, unsigned *npanel);
EXTERNC void cal_dskpfa_(double   *Pfa, char *uplo, unsigned *n, double   *A, unsigned *ldA, unsigned *inv,
        double   *Sp1, double   *Sp2, double   *Sp3, double   *Sp4, double   *Sp5, unsigned *npanel);
EXTERNC void cal_cskpfa_(scomplex *Pfa, char *uplo, unsigned *n, scomplex *A, unsigned *ldA, unsigned *inv,
        scomplex *Sp1, scomplex *Sp2, scomplex *Sp3, scomplex *Sp4, scomplex *Sp5, unsigned *npanel);
EXTERNC void cal_zskpfa_(dcomplex *Pfa, char *uplo, unsigned *n, dcomplex *A, unsigned *ldA, unsigned *inv,
        dcomplex *Sp1, dcomplex *Sp2, dcomplex *Sp3, dcomplex *Sp4, dcomplex *Sp5, unsigned *npanel);

