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
typedef std::complex<float>  ccscmplx;
typedef std::complex<double> ccdcmplx;
#else
#define EXTERNC
#include <complex.h>
typedef float  complex ccscmplx;
typedef double complex ccdcmplx;
#endif

// Simplified interface.
EXTERNC void sskpfa(float    *Pfa, char uplo, unsigned n, float    *A, unsigned ldA, unsigned inv);
EXTERNC void dskpfa(double   *Pfa, char uplo, unsigned n, double   *A, unsigned ldA, unsigned inv);
EXTERNC void cskpfa(ccscmplx *Pfa, char uplo, unsigned n, ccscmplx *A, unsigned ldA, unsigned inv);
EXTERNC void zskpfa(ccdcmplx *Pfa, char uplo, unsigned n, ccdcmplx *A, unsigned ldA, unsigned inv);

// Full interface, callable from both C and Fortran.
EXTERNC void cal_sskpfa_(float    *Pfa, char *uplo, unsigned *n, float    *A, unsigned *ldA, unsigned *inv,
        float    *Sp1, float    *Sp2, float    *Sp3, float    *Sp4, float    *Sp5, unsigned *npanel);
EXTERNC void cal_dskpfa_(double   *Pfa, char *uplo, unsigned *n, double   *A, unsigned *ldA, unsigned *inv,
        double   *Sp1, double   *Sp2, double   *Sp3, double   *Sp4, double   *Sp5, unsigned *npanel);
EXTERNC void cal_cskpfa_(ccscmplx *Pfa, char *uplo, unsigned *n, ccscmplx *A, unsigned *ldA, unsigned *inv,
        ccscmplx *Sp1, ccscmplx *Sp2, ccscmplx *Sp3, ccscmplx *Sp4, ccscmplx *Sp5, unsigned *npanel);
EXTERNC void cal_zskpfa_(ccdcmplx *Pfa, char *uplo, unsigned *n, ccdcmplx *A, unsigned *ldA, unsigned *inv,
        ccdcmplx *Sp1, ccdcmplx *Sp2, ccdcmplx *Sp3, ccdcmplx *Sp4, ccdcmplx *Sp5, unsigned *npanel);

