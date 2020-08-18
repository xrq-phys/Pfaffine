/**
 * \file blalink.hh
 * C++ to C type wrapper for BLIS type-apis.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#pragma once
#include <complex>
typedef std::complex<float>  ccscmplx;
typedef std::complex<double> ccdcmplx;
// BLIS definitions.
#include "blis.h"
// BLAS interface.
#include "blalink_fort.h"


// gemm
template <typename T>
inline void gemm(trans_t transa, trans_t transb,
                 dim_t m, dim_t n, dim_t k,
                 T alpha,
                 T *a, inc_t lda,
                 T *b, inc_t ldb,
                 T beta,
                 T *c, inc_t ldc);
template <> inline void gemm<float>(trans_t transa, trans_t transb, dim_t m, dim_t n, dim_t k, float alpha, float *a, inc_t lda, float *b, inc_t ldb, float beta, float *c, inc_t ldc)
{ char ta = trans2char(transa), tb = trans2char(transb); sgemm_(&ta, &tb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); }
template <> inline void gemm<double>(trans_t transa, trans_t transb, dim_t m, dim_t n, dim_t k, double alpha, double *a, inc_t lda, double *b, inc_t ldb, double beta, double *c, inc_t ldc)
{ char ta = trans2char(transa), tb = trans2char(transb); dgemm_(&ta, &tb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); }
template <> inline void gemm<ccscmplx>(trans_t transa, trans_t transb, dim_t m, dim_t n, dim_t k, ccscmplx alpha, ccscmplx *a, inc_t lda, ccscmplx *b, inc_t ldb, ccscmplx beta, ccscmplx *c, inc_t ldc)
{ char ta = trans2char(transa), tb = trans2char(transb); cgemm_(&ta, &tb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); }
template <> inline void gemm<ccdcmplx>(trans_t transa, trans_t transb, dim_t m, dim_t n, dim_t k, ccdcmplx alpha, ccdcmplx *a, inc_t lda, ccdcmplx *b, inc_t ldb, ccdcmplx beta, ccdcmplx *c, inc_t ldc)
{ char ta = trans2char(transa), tb = trans2char(transb); zgemm_(&ta, &tb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); }


// ger
template <typename T>
inline void ger(dim_t m, dim_t n,
                T alpha,
                T *x, inc_t incx,
                T *y, inc_t incy,
                T *a, inc_t lda);
template <> inline void ger<float>(dim_t m, dim_t n, float alpha, float *x, inc_t incx, float *y, inc_t incy, float *a, inc_t lda)
{ sger_(&m, &n, &alpha, x, &incx, y, &incy, a, &lda); }
template <> inline void ger<double>(dim_t m, dim_t n, double alpha, double *x, inc_t incx, double *y, inc_t incy, double *a, inc_t lda)
{ dger_(&m, &n, &alpha, x, &incx, y, &incy, a, &lda); }
template <> inline void ger<ccscmplx>(dim_t m, dim_t n, ccscmplx alpha, ccscmplx *x, inc_t incx, ccscmplx *y, inc_t incy, ccscmplx *a, inc_t lda)
{ cgeru_(&m, &n, &alpha, x, &incx, y, &incy, a, &lda); }
template <> inline void ger<ccdcmplx>(dim_t m, dim_t n, ccdcmplx alpha, ccdcmplx *x, inc_t incx, ccdcmplx *y, inc_t incy, ccdcmplx *a, inc_t lda)
{ zgeru_(&m, &n, &alpha, x, &incx, y, &incy, a, &lda); }


// gemv
template <typename T>
inline void gemv(trans_t trans,
                 dim_t m, dim_t n,
                 T alpha,
                 T *a, inc_t lda,
                 T *x, inc_t incx,
                 T beta,
                 T *y, inc_t incy);
template <> inline void gemv<float>(trans_t trans, dim_t m, dim_t n, float alpha, float *a, inc_t lda, float *x, inc_t incx, float beta, float *y, inc_t incy)
{ char t = trans2char(trans); sgemv_(&t, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy); }
template <> inline void gemv<double>(trans_t trans, dim_t m, dim_t n, double alpha, double *a, inc_t lda, double *x, inc_t incx, double beta, double *y, inc_t incy)
{ char t = trans2char(trans); dgemv_(&t, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy); }
template <> inline void gemv<ccscmplx>(trans_t trans, dim_t m, dim_t n, ccscmplx alpha, ccscmplx *a, inc_t lda, ccscmplx *x, inc_t incx, ccscmplx beta, ccscmplx *y, inc_t incy)
{ char t = trans2char(trans); cgemv_(&t, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy); }
template <> inline void gemv<ccdcmplx>(trans_t trans, dim_t m, dim_t n, ccdcmplx alpha, ccdcmplx *a, inc_t lda, ccdcmplx *x, inc_t incx, ccdcmplx beta, ccdcmplx *y, inc_t incy)
{ char t = trans2char(trans); zgemv_(&t, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy); }


// swap
template <typename T>
inline void swap(dim_t n,
                 T *x, inc_t incx,
                 T *y, inc_t incy);
template <> inline void swap<float>(dim_t n, float *x, inc_t incx, float *y, inc_t incy)
{ sswap_(&n, x, &incx, y, &incy); }
template <> inline void swap<double>(dim_t n, double *x, inc_t incx, double *y, inc_t incy)
{ dswap_(&n, x, &incx, y, &incy); }
template <> inline void swap<ccscmplx>(dim_t n, ccscmplx *x, inc_t incx, ccscmplx *y, inc_t incy)
{ cswap_(&n, x, &incx, y, &incy); }
template <> inline void swap<ccdcmplx>(dim_t n, ccdcmplx *x, inc_t incx, ccdcmplx *y, inc_t incy)
{ zswap_(&n, x, &incx, y, &incy); }

// axpy
template <typename T>
inline void axpy(dim_t n,
                 T alpha,
                 T *x, inc_t incx,
                 T *y, inc_t incy);
template <> inline void axpy<float>(dim_t n, float alpha, float *x, inc_t incx, float *y, inc_t incy)
{ saxpy_(&n, &alpha, x, &incx, y, &incy); }
template <> inline void axpy<double>(dim_t n, double alpha, double *x, inc_t incx, double *y, inc_t incy)
{ daxpy_(&n, &alpha, x, &incx, y, &incy); }
template <> inline void axpy<ccscmplx>(dim_t n, ccscmplx alpha, ccscmplx *x, inc_t incx, ccscmplx *y, inc_t incy)
{ caxpy_(&n, &alpha, x, &incx, y, &incy); }
template <> inline void axpy<ccdcmplx>(dim_t n, ccdcmplx alpha, ccdcmplx *x, inc_t incx, ccdcmplx *y, inc_t incy)
{ zaxpy_(&n, &alpha, x, &incx, y, &incy); }

// dot
template <typename T>
inline T dot(dim_t n,
             T *sx, inc_t incx,
             T *sy, inc_t incy);
template <> inline float dot<float>(dim_t n, float *sx, inc_t incx, float *sy, inc_t incy)
{ return sdot_(&n, sx, &incx, sy, &incy); }
template <> inline double dot<double>(dim_t n, double *sx, inc_t incx, double *sy, inc_t incy)
{ return ddot_(&n, sx, &incx, sy, &incy); }
template <> inline ccscmplx dot<ccscmplx>(dim_t n, ccscmplx *sx, inc_t incx, ccscmplx *sy, inc_t incy)
{ scomplex rho = cdotc_(&n, sx, &incx, sy, &incy); return *((ccscmplx *)&rho); }
template <> inline ccdcmplx dot<ccdcmplx>(dim_t n, ccdcmplx *sx, inc_t incx, ccdcmplx *sy, inc_t incy)
{ dcomplex rho = zdotc_(&n, sx, &incx, sy, &incy); return *((ccdcmplx *)&rho); }


// [extension] skr2k
template <typename T>
inline void skr2k(uplo_t uploc,
                  trans_t transab,
                  dim_t m, dim_t k,
                  T alpha,
                  T *a, inc_t lda,
                  T *b, inc_t ldb,
                  T beta,
                  T *c, inc_t ldc);
template <> inline void skr2k<float>(uplo_t uploc, trans_t transab, dim_t m, dim_t k, float alpha, float *a, inc_t lda, float *b, inc_t ldb, float beta, float *c, inc_t ldc)
{ bli_sskr2k(uploc, transab, transab, m, k, &alpha, a, 1, lda, b, 1, ldb, &beta, c, 1, ldc); }
template <> inline void skr2k<double>(uplo_t uploc, trans_t transab, dim_t m, dim_t k, double alpha, double *a, inc_t lda, double *b, inc_t ldb, double beta, double *c, inc_t ldc)
{ bli_dskr2k(uploc, transab, transab, m, k, &alpha, a, 1, lda, b, 1, ldb, &beta, c, 1, ldc); }
template <> inline void skr2k<ccscmplx>(uplo_t uploc, trans_t transab, dim_t m, dim_t k, ccscmplx alpha, ccscmplx *a, inc_t lda, ccscmplx *b, inc_t ldb, ccscmplx beta, ccscmplx *c, inc_t ldc)
{ bli_cskr2k(uploc, transab, transab, m, k, (scomplex *)&alpha, (scomplex *)a, 1, lda, (scomplex *)b, 1, ldb, (scomplex *)&beta, (scomplex *)c, 1, ldc); }
template <> inline void skr2k<ccdcmplx>(uplo_t uploc, trans_t transab, dim_t m, dim_t k, ccdcmplx alpha, ccdcmplx *a, inc_t lda, ccdcmplx *b, inc_t ldb, ccdcmplx beta, ccdcmplx *c, inc_t ldc)
{ bli_zskr2k(uploc, transab, transab, m, k, (dcomplex *)&alpha, (dcomplex *)a, 1, lda, (dcomplex *)b, 1, ldb, (dcomplex *)&beta, (dcomplex *)c, 1, ldc); }
