/**
 * \file kersel.hh
 * gemm kernel selector for skr2k.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#pragma once
#include <complex>
typedef std::complex<float>  scomplex;
typedef std::complex<double> dcomplex;

// if has optimized microkernel available.
template <typename T> inline unsigned mker_available(void);
template <> inline unsigned mker_available<float>   (void) { return 0; }
template <> inline unsigned mker_available<double>  (void)
#if defined(_Sandy)
{ return 1; }
#else
{ return 0; }
#endif
template <> inline unsigned mker_available<scomplex>(void) { return 0; }
template <> inline unsigned mker_available<dcomplex>(void) { return 0; }
// block size
template <typename T> inline void set_blk_size(unsigned *mr, unsigned *nr);
template <> inline void set_blk_size<float>   (unsigned *mr, unsigned *nr) { *mr = 6; *nr = 8; }
template <> inline void set_blk_size<double>  (unsigned *mr, unsigned *nr) { *mr = 6; *nr = 8; }
template <> inline void set_blk_size<scomplex>(unsigned *mr, unsigned *nr) { *mr = 4; *nr = 4; }
template <> inline void set_blk_size<dcomplex>(unsigned *mr, unsigned *nr) { *mr = 4; *nr = 4; }

// defined kernels.
#if defined(_Sandy)
extern "C" void udgemmn(unsigned k, double *alpha, double *pakA, double *pakB, double *beta, double *C, unsigned ldC);
extern "C" void udgemmt(unsigned k, double *alpha, double *pakA, double *pakB, double *beta, double *C, unsigned ldC);
#endif

// gemm selectors.
// mr * nr kernels.
inline void ugemmn(unsigned k, float *alpha, float *pakA, float *pakB, float *beta, float *C, unsigned ldC)
{ std::_Exit(EXIT_FAILURE); }
inline void ugemmn(unsigned k, double *alpha, double *pakA, double *pakB, double *beta, double *C, unsigned ldC)
#if defined(_Sandy)
{ udgemmn(k, alpha, pakA, pakB, beta, C, ldC); }
#else
{ std::_Exit(EXIT_FAILURE); }
#endif
inline void ugemmn(unsigned k, scomplex *alpha, scomplex *pakA, scomplex *pakB, scomplex *beta, scomplex *C, unsigned ldC)
{ std::_Exit(EXIT_FAILURE); }
inline void ugemmn(unsigned k, dcomplex *alpha, dcomplex *pakA, dcomplex *pakB, dcomplex *beta, dcomplex *C, unsigned ldC)
{ std::_Exit(EXIT_FAILURE); }
// nr * mr kernels.
inline void ugemmt(unsigned k, float *alpha, float *pakA, float *pakB, float *beta, float *C, unsigned ldC)
{ std::_Exit(EXIT_FAILURE); }
inline void ugemmt(unsigned k, double *alpha, double *pakA, double *pakB, double *beta, double *C, unsigned ldC)
#if defined(_Sandy)
{ udgemmn(k, alpha, pakA, pakB, beta, C, ldC); /*Cases udgemmt is the same as udgemmn.*/ }
#else
{ std::_Exit(EXIT_FAILURE); }
#endif
inline void ugemmt(unsigned k, scomplex *alpha, scomplex *pakA, scomplex *pakB, scomplex *beta, scomplex *C, unsigned ldC)
{ std::_Exit(EXIT_FAILURE); }
inline void ugemmt(unsigned k, dcomplex *alpha, dcomplex *pakA, dcomplex *pakB, dcomplex *beta, dcomplex *C, unsigned ldC)
{ std::_Exit(EXIT_FAILURE); }

