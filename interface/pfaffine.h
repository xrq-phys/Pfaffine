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

#ifdef EXPAND_NAME
#undef EXPAND_NAME
#endif
#define EXPAND_NAME( typechar, funcname ) typechar##funcname##_

// Simplified interface.
#ifdef PASTE_DEF
#undef PASTE_DEF
#endif
#define PASTE_DEF( typename, typechar ) \
    EXTERNC void EXPAND_NAME( typechar, skpfa ) \
        (char uplo, unsigned n, float *A, signed ldA, unsigned inv, typename *dPfa, signed *info);

PASTE_DEF( float,    s )
PASTE_DEF( double,   d )
PASTE_DEF( ccscmplx, c )
PASTE_DEF( ccdcmplx, z )
#undef EXPAND_NAME
#undef PASTE_DEF

// Full interface, callable from both C and Fortran.
#define EXPAND_NAME( typechar, funcname ) cal_##typechar##funcname##_
#define PASTE_DEF( typename, typechar ) \
    EXTERNC void EXPAND_NAME( typechar, skpfa ) \
        (char *uplo,  \
         unsigned *n, \
         typename *_A, signed *ldA, \
         typename *_G, signed *ldG, \
         signed *iPov,   \
         signed *inv,    \
         typename *dPfa, signed *info, \
         typename *_Work, unsigned *lWork);

PASTE_DEF( float,    s )
PASTE_DEF( double,   d )
PASTE_DEF( ccscmplx, c )
PASTE_DEF( ccdcmplx, z )
#undef EXPAND_NAME
#undef PASTE_DEF
