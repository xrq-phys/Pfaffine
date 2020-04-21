/**
 * \file kersel.hh
 * gemm custom-sized kernel selector for skr2k.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#pragma once
#include <complex>
typedef std::complex<float>  scomplex;
typedef std::complex<double> dcomplex;

// availability.
template <typename T> inline unsigned extker_available(void);
template <> inline unsigned extker_available<float>   (void)
{ return 0; }
template <> inline unsigned extker_available<double>  (void)
#if defined(_SVE)
{ return 1; }
#else
{ return 0; }
#endif
template <> inline unsigned extker_available<scomplex>(void)
{ return 0; }
template <> inline unsigned extker_available<dcomplex>(void)
{ return 0; }

extern "C" {
#if defined(_SVE)
unsigned udgemmext(unsigned m, unsigned n, unsigned k,
                   double *Alpha_, double *A, unsigned ldA, double *B, unsigned ldB,
                   double *Beta_, double *C, unsigned ldC);
#endif
}

unsigned ugemmext(unsigned m, unsigned n, unsigned k,
                  float *Alpha_, 
                  float *A, unsigned ldA, 
                  float *B, unsigned ldB,
                  float *Beta_, 
                  float *C, unsigned ldC)
{ const static unsigned info_err = 1; return info_err; }
unsigned ugemmext(unsigned m, unsigned n, unsigned k,
                  double *Alpha_, 
                  double *A, unsigned ldA, 
                  double *B, unsigned ldB,
                  double *Beta_, 
                  double *C, unsigned ldC)
#if defined(_SVE)
{ return udgemmext(m, n, k, Alpha_, A, ldA, B, ldB, Beta_, C, ldC); }
#else
{ const static unsigned info_err = 1; return info_err; }
#endif
unsigned ugemmext(unsigned m, unsigned n, unsigned k,
                  scomplex *Alpha_, 
                  scomplex *A, unsigned ldA, 
                  scomplex *B, unsigned ldB,
                  scomplex *Beta_, 
                  scomplex *C, unsigned ldC)
{ const static unsigned info_err = 1; return info_err; }
unsigned ugemmext(unsigned m, unsigned n, unsigned k,
                  dcomplex *Alpha_, 
                  dcomplex *A, unsigned ldA, 
                  dcomplex *B, unsigned ldB,
                  dcomplex *Beta_, 
                  dcomplex *C, unsigned ldC)
{ const static unsigned info_err = 1; return info_err; }

