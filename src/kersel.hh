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

// Memory alignment block size.
// TODO: Determine for SVE by vector bits.
#if defined(_SkylakeX) || defined(_SVE)
const unsigned align_blk = 64;
#else
const unsigned align_blk = 32;
#endif

#if defined(_SVE)
// Vector length query for SVE.
extern "C" unsigned dvecln_iso(void);
extern "C" unsigned zgemm_sve_mr(void);
extern "C" unsigned zgemm_sve_nr(void);
#endif

// if has optimized microkernel available.
template <typename T> inline unsigned mker_available(void);
template <> inline unsigned mker_available<float>   (void)
#if defined(_Neon) || defined(_Haswell)
{ return 1; }
#else
{ return 0; }
#endif
template <> inline unsigned mker_available<double>  (void)
#if defined(_Sandy) || defined(_Neon) || defined(_Haswell) || defined(_SVE)
{ return 1; }
#else
{ return 0; }
#endif
template <> inline unsigned mker_available<scomplex>(void)
#if defined(_Haswell)
{ return 1; }
#else
{ return 0; }
#endif
template <> inline unsigned mker_available<dcomplex>(void)
#if defined(_Haswell) || defined(_SVE)
{ return 1; }
#else
{ return 0; }
#endif
// if a complex type has corresponding real kernel.
// roadmap: utillize real kernels to compute complex.
template <typename T> inline unsigned elemker_available(void);
template <> inline unsigned elemker_available<scomplex>(void) { return mker_available<float> (); }
template <> inline unsigned elemker_available<dcomplex>(void) { return mker_available<double>(); }
template <> inline unsigned elemker_available<float >(void) { return 0; }
template <> inline unsigned elemker_available<double>(void) { return 0; }
// block size
template <typename T> inline void set_blk_size(unsigned *mr, unsigned *nr);
template <> inline void set_blk_size<float>   (unsigned *mr, unsigned *nr)
#if defined(_Haswell)
{ *mr = 6; *nr = 16; }
#else
{ *mr = 8; *nr = 12; }
#endif
template <> inline void set_blk_size<double>  (unsigned *mr, unsigned *nr)
#if defined(_Sandy)
{ *mr = 4; *nr = 4; }
#elif defined(_SVE)
{ *mr = dvecln_iso(); *nr = dvecln_iso(); }
#elif defined(_Haswell)
#if defined(_SkylakeX)
{ *mr = 14; *nr = 16; }
#else
{ *mr = 6; *nr = 8; }
#endif
#else
{ *mr = 8; *nr = 8; }
#endif
template <> inline void set_blk_size<scomplex>(unsigned *mr, unsigned *nr) { *mr = 3; *nr = 8; }
template <> inline void set_blk_size<dcomplex>(unsigned *mr, unsigned *nr)
#if defined(_SVE)
{ *mr = zgemm_sve_mr(); *nr = zgemm_sve_nr(); }
#else
{ *mr = 3; *nr = 4; }
#endif

// defined kernels.
extern "C" {
#if defined(_Haswell) || defined(_Neon)
void usgemmn(unsigned k, float *alpha, float *pakA, float *pakB, float *beta, float *C, unsigned ldC);
void usgemmt(unsigned k, float *alpha, float *pakA, float *pakB, float *beta, float *C, unsigned ldC);
#endif
#if defined(_Haswell) || defined(_Sandy) || defined(_Neon) || defined(_SVE)
void udgemmn(unsigned k, double *alpha, double *pakA, double *pakB, double *beta, double *C, unsigned ldC);
void udgemmt(unsigned k, double *alpha, double *pakA, double *pakB, double *beta, double *C, unsigned ldC);
#endif
#if defined(_Haswell)
void ucgemmn(unsigned k, scomplex *alpha, scomplex *pakA, scomplex *pakB, scomplex *beta, scomplex *C, unsigned ldC);
void ucgemmt(unsigned k, scomplex *alpha, scomplex *pakA, scomplex *pakB, scomplex *beta, scomplex *C, unsigned ldC);
#endif
#if defined(_Haswell) || defined(_SVE)
void uzgemmn(unsigned k, dcomplex *alpha, dcomplex *pakA, dcomplex *pakB, dcomplex *beta, dcomplex *C, unsigned ldC);
void uzgemmt(unsigned k, dcomplex *alpha, dcomplex *pakA, dcomplex *pakB, dcomplex *beta, dcomplex *C, unsigned ldC);
#endif
}

// gemm selectors.
// mr * nr kernels.
inline void ugemmn(unsigned k, float *alpha, float *pakA, float *pakB, float *beta, float *C, unsigned ldC)
#if defined(_Haswell) || defined(_Neon)
{ usgemmn(k, alpha, pakA, pakB, beta, C, ldC); }
#else
{ std::_Exit(EXIT_FAILURE); }
#endif
inline void ugemmn(unsigned k, double *alpha, double *pakA, double *pakB, double *beta, double *C, unsigned ldC)
#if defined(_Haswell) || defined(_Sandy) || defined(_Neon) || defined(_SVE)
{ udgemmn(k, alpha, pakA, pakB, beta, C, ldC); }
#else
{ std::_Exit(EXIT_FAILURE); }
#endif
inline void ugemmn(unsigned k, scomplex *alpha, scomplex *pakA, scomplex *pakB, scomplex *beta, scomplex *C, unsigned ldC)
#if defined(_Haswell)
{ ucgemmn(k, alpha, pakA, pakB, beta, C, ldC); }
#else
{ std::_Exit(EXIT_FAILURE); }
#endif
inline void ugemmn(unsigned k, dcomplex *alpha, dcomplex *pakA, dcomplex *pakB, dcomplex *beta, dcomplex *C, unsigned ldC)
#if defined(_Haswell) || defined(_SVE)
{ uzgemmn(k, alpha, pakA, pakB, beta, C, ldC); }
#else
{ std::_Exit(EXIT_FAILURE); }
#endif
// nr * mr kernels.
inline void ugemmt(unsigned k, float *alpha, float *pakA, float *pakB, float *beta, float *C, unsigned ldC)
#if defined(_Haswell) || defined(_Neon)
{ usgemmt(k, alpha, pakA, pakB, beta, C, ldC); }
#else
{ std::_Exit(EXIT_FAILURE); }
#endif
inline void ugemmt(unsigned k, double *alpha, double *pakA, double *pakB, double *beta, double *C, unsigned ldC)
#if defined(_Haswell) || defined(_Sandy) || defined(_Neon) || defined(_SVE)
{ udgemmt(k, alpha, pakA, pakB, beta, C, ldC); }
#else
{ std::_Exit(EXIT_FAILURE); }
#endif
inline void ugemmt(unsigned k, scomplex *alpha, scomplex *pakA, scomplex *pakB, scomplex *beta, scomplex *C, unsigned ldC)
#if defined(_Haswell)
{ ucgemmt(k, alpha, pakA, pakB, beta, C, ldC); }
#else
{ std::_Exit(EXIT_FAILURE); }
#endif
inline void ugemmt(unsigned k, dcomplex *alpha, dcomplex *pakA, dcomplex *pakB, dcomplex *beta, dcomplex *C, unsigned ldC)
#if defined(_Haswell) || defined(_SVE)
{ uzgemmt(k, alpha, pakA, pakB, beta, C, ldC); }
#else
{ std::_Exit(EXIT_FAILURE); }
#endif

