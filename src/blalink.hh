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

// error info
typedef enum {
    Pfaffine_SIGN_ERR = 1,
    Pfaffine_OUT_OF_BOUND,
    Pfaffine_BAD_SCRATCHPAD,
    Pfaffine_BAD_REPRESENTATION,
    Pfaffine_INTEGER_OVERFLOW,
    Pfaffine_DOUBLE_NAN_DETECTED,
    Pfaffine_NOT_IMPLEMNTED,
    Pfaffine_NUM_ERROR_TYPE
} skpfa_error_t;
inline signed err_info(signed type, signed pos)
{ return type * Pfaffine_NUM_ERROR_TYPE + pos; }

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
{ bli_sgemm(transa, transb, m, n, k, &alpha, a, 1, lda, b, 1, ldb, &beta, c, 1, ldc); }
template <> inline void gemm<double>(trans_t transa, trans_t transb, dim_t m, dim_t n, dim_t k, double alpha, double *a, inc_t lda, double *b, inc_t ldb, double beta, double *c, inc_t ldc)
{ bli_dgemm(transa, transb, m, n, k, &alpha, a, 1, lda, b, 1, ldb, &beta, c, 1, ldc); }
template <> inline void gemm<ccscmplx>(trans_t transa, trans_t transb, dim_t m, dim_t n, dim_t k, ccscmplx alpha, ccscmplx *a, inc_t lda, ccscmplx *b, inc_t ldb, ccscmplx beta, ccscmplx *c, inc_t ldc)
{ bli_cgemm(transa, transb, m, n, k, (scomplex *)&alpha, (scomplex *)a, 1, lda, (scomplex *)b, 1, ldb, (scomplex *)&beta, (scomplex *)c, 1, ldc); }
template <> inline void gemm<ccdcmplx>(trans_t transa, trans_t transb, dim_t m, dim_t n, dim_t k, ccdcmplx alpha, ccdcmplx *a, inc_t lda, ccdcmplx *b, inc_t ldb, ccdcmplx beta, ccdcmplx *c, inc_t ldc)
{ bli_zgemm(transa, transb, m, n, k, (dcomplex *)&alpha, (dcomplex *)a, 1, lda, (dcomplex *)b, 1, ldb, (dcomplex *)&beta, (dcomplex *)c, 1, ldc); }


// ger
template <typename T>
inline void ger(dim_t m, dim_t n,
                T alpha,
                T *x, inc_t incx,
                T *y, inc_t incy,
                T *a, inc_t lda);
template <> inline void ger<float>(dim_t m, dim_t n, float alpha, float *x, inc_t incx, float *y, inc_t incy, float *a, inc_t lda)
{ bli_sger(BLIS_NO_CONJUGATE, BLIS_NO_CONJUGATE, m, n, &alpha, x, incx, y, incy, a, 1, lda); }
template <> inline void ger<double>(dim_t m, dim_t n, double alpha, double *x, inc_t incx, double *y, inc_t incy, double *a, inc_t lda)
{ bli_dger(BLIS_NO_CONJUGATE, BLIS_NO_CONJUGATE, m, n, &alpha, x, incx, y, incy, a, 1, lda); }
template <> inline void ger<ccscmplx>(dim_t m, dim_t n, ccscmplx alpha, ccscmplx *x, inc_t incx, ccscmplx *y, inc_t incy, ccscmplx *a, inc_t lda)
{ bli_cger(BLIS_NO_CONJUGATE, BLIS_NO_CONJUGATE, m, n, (scomplex *)&alpha, (scomplex *)x, incx, (scomplex *)y, incy, (scomplex *)a, 1, lda); }
template <> inline void ger<ccdcmplx>(dim_t m, dim_t n, ccdcmplx alpha, ccdcmplx *x, inc_t incx, ccdcmplx *y, inc_t incy, ccdcmplx *a, inc_t lda)
{ bli_zger(BLIS_NO_CONJUGATE, BLIS_NO_CONJUGATE, m, n, (dcomplex *)&alpha, (dcomplex *)x, incx, (dcomplex *)y, incy, (dcomplex *)a, 1, lda); }


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
{ bli_sgemv(trans, BLIS_NO_CONJUGATE, m, n, &alpha, a, 1, lda, x, incx, &beta, y, incy); }
template <> inline void gemv<double>(trans_t trans, dim_t m, dim_t n, double alpha, double *a, inc_t lda, double *x, inc_t incx, double beta, double *y, inc_t incy)
{ bli_dgemv(trans, BLIS_NO_CONJUGATE, m, n, &alpha, a, 1, lda, x, incx, &beta, y, incy); }
template <> inline void gemv<ccscmplx>(trans_t trans, dim_t m, dim_t n, ccscmplx alpha, ccscmplx *a, inc_t lda, ccscmplx *x, inc_t incx, ccscmplx beta, ccscmplx *y, inc_t incy)
{ bli_cgemv(trans, BLIS_NO_CONJUGATE, m, n, (scomplex *)&alpha, (scomplex *)a, 1, lda, (scomplex *)x, incx, (scomplex *)&beta, (scomplex *)y, incy); }
template <> inline void gemv<ccdcmplx>(trans_t trans, dim_t m, dim_t n, ccdcmplx alpha, ccdcmplx *a, inc_t lda, ccdcmplx *x, inc_t incx, ccdcmplx beta, ccdcmplx *y, inc_t incy)
{ bli_zgemv(trans, BLIS_NO_CONJUGATE, m, n, (dcomplex *)&alpha, (dcomplex *)a, 1, lda, (dcomplex *)x, incx, (dcomplex *)&beta, (dcomplex *)y, incy); }


// trmm
template <typename T>
inline void trmm(side_t sidea, 
                 uplo_t uploa, 
                 trans_t transa, 
                 dim_t m, dim_t n,
                 T alpha, 
                 T *a, inc_t lda, 
                 T *b, inc_t ldb);
template <> inline void trmm<float>(side_t sidea, uplo_t uploa, trans_t transa, dim_t m, dim_t n, float alpha, float *a, inc_t lda, float *b, inc_t ldb) 
{ bli_strmm(sidea, uploa, transa, BLIS_NONUNIT_DIAG, m, n, &alpha, a, 1, lda, b, 1, ldb); }
template <> inline void trmm<double>(side_t sidea, uplo_t uploa, trans_t transa, dim_t m, dim_t n, double alpha, double *a, inc_t lda, double *b, inc_t ldb) 
{ bli_dtrmm(sidea, uploa, transa, BLIS_NONUNIT_DIAG, m, n, &alpha, a, 1, lda, b, 1, ldb); }
template <> inline void trmm<ccscmplx>(side_t sidea, uplo_t uploa, trans_t transa, dim_t m, dim_t n, ccscmplx alpha, ccscmplx *a, inc_t lda, ccscmplx *b, inc_t ldb) 
{ bli_ctrmm(sidea, uploa, transa, BLIS_NONUNIT_DIAG, m, n, (scomplex *)&alpha, (scomplex *)a, 1, lda, (scomplex *)b, 1, ldb); }
template <> inline void trmm<ccdcmplx>(side_t sidea, uplo_t uploa, trans_t transa, dim_t m, dim_t n, ccdcmplx alpha, ccdcmplx *a, inc_t lda, ccdcmplx *b, inc_t ldb) 
{ bli_ztrmm(sidea, uploa, transa, BLIS_NONUNIT_DIAG, m, n, (dcomplex *)&alpha, (dcomplex *)a, 1, lda, (dcomplex *)b, 1, ldb); }

// trmv
template <typename T>
inline void trmv(uplo_t uploa,
                 trans_t transa,
                 dim_t m,
                 T alpha,
                 T *a, inc_t lda,
                 T *x, inc_t incx);
template <> inline void trmv<float>(uplo_t uploa, trans_t transa, dim_t m, float alpha, float *a, inc_t lda, float *x, inc_t incx)
{ bli_strmv(uploa, transa, BLIS_NONUNIT_DIAG, m, &alpha, a, 1, lda, x, incx); }
template <> inline void trmv<double>(uplo_t uploa, trans_t transa, dim_t m, double alpha, double *a, inc_t lda, double *x, inc_t incx)
{ bli_dtrmv(uploa, transa, BLIS_NONUNIT_DIAG, m, &alpha, a, 1, lda, x, incx); }
template <> inline void trmv<ccscmplx>(uplo_t uploa, trans_t transa, dim_t m, ccscmplx alpha, ccscmplx *a, inc_t lda, ccscmplx *x, inc_t incx)
{ bli_ctrmv(uploa, transa, BLIS_NONUNIT_DIAG, m, (scomplex *)&alpha, (scomplex *)a, 1, lda, (scomplex *)x, incx); }
template <> inline void trmv<ccdcmplx>(uplo_t uploa, trans_t transa, dim_t m, ccdcmplx alpha, ccdcmplx *a, inc_t lda, ccdcmplx *x, inc_t incx)
{ bli_ztrmv(uploa, transa, BLIS_NONUNIT_DIAG, m, (dcomplex *)&alpha, (dcomplex *)a, 1, lda, (dcomplex *)x, incx); }

// swap
template <typename T>
inline void swap(dim_t n,
                 T *x, inc_t incx,
                 T *y, inc_t incy);
template <> inline void swap<float>(dim_t n, float *x, inc_t incx, float *y, inc_t incy)
{ bli_sswapv(n, x, incx, y, incy); }
template <> inline void swap<double>(dim_t n, double *x, inc_t incx, double *y, inc_t incy)
{ bli_dswapv(n, x, incx, y, incy); }
template <> inline void swap<ccscmplx>(dim_t n, ccscmplx *x, inc_t incx, ccscmplx *y, inc_t incy)
{ bli_cswapv(n, (scomplex *)x, incx, (scomplex *)y, incy); }
template <> inline void swap<ccdcmplx>(dim_t n, ccdcmplx *x, inc_t incx, ccdcmplx *y, inc_t incy)
{ bli_zswapv(n, (dcomplex *)x, incx, (dcomplex *)y, incy); }

// axpy
template <typename T>
inline void axpy(dim_t n,
                 T alpha,
                 T *x, inc_t incx,
                 T *y, inc_t incy);
template <> inline void axpy<float>(dim_t n, float alpha, float *x, inc_t incx, float *y, inc_t incy)
{ bli_saxpyv(BLIS_NO_CONJUGATE, n, &alpha, x, incx, y, incy); }
template <> inline void axpy<double>(dim_t n, double alpha, double *x, inc_t incx, double *y, inc_t incy)
{ bli_daxpyv(BLIS_NO_CONJUGATE, n, &alpha, x, incx, y, incy); }
template <> inline void axpy<ccscmplx>(dim_t n, ccscmplx alpha, ccscmplx *x, inc_t incx, ccscmplx *y, inc_t incy)
{ bli_caxpyv(BLIS_NO_CONJUGATE, n, (scomplex *)&alpha, (scomplex *)x, incx, (scomplex *)y, incy); }
template <> inline void axpy<ccdcmplx>(dim_t n, ccdcmplx alpha, ccdcmplx *x, inc_t incx, ccdcmplx *y, inc_t incy)
{ bli_zaxpyv(BLIS_NO_CONJUGATE, n, (dcomplex *)&alpha, (dcomplex *)x, incx, (dcomplex *)y, incy); }

// dot
template <typename T>
inline T dot(dim_t n,
             T *sx, inc_t incx,
             T *sy, inc_t incy);
template <> inline float dot<float>(dim_t n, float *sx, inc_t incx, float *sy, inc_t incy)
{ float rho; bli_sdotv(BLIS_NO_CONJUGATE, BLIS_NO_CONJUGATE, n, sx, incx, sy, incy, &rho); return rho; }
template <> inline double dot<double>(dim_t n, double *sx, inc_t incx, double *sy, inc_t incy)
{ double rho; bli_ddotv(BLIS_NO_CONJUGATE, BLIS_NO_CONJUGATE, n, sx, incx, sy, incy, &rho); return rho; }
template <> inline ccscmplx dot<ccscmplx>(dim_t n, ccscmplx *sx, inc_t incx, ccscmplx *sy, inc_t incy)
{ ccscmplx rho; bli_cdotv(BLIS_NO_CONJUGATE, BLIS_NO_CONJUGATE, n, (scomplex *)sx, incx, (scomplex *)sy, incy, (scomplex *)&rho); return rho; }
template <> inline ccdcmplx dot<ccdcmplx>(dim_t n, ccdcmplx *sx, inc_t incx, ccdcmplx *sy, inc_t incy)
{ ccdcmplx rho; bli_zdotv(BLIS_NO_CONJUGATE, BLIS_NO_CONJUGATE, n, (dcomplex *)sx, incx, (dcomplex *)sy, incy, (dcomplex *)&rho); return rho; }


// [BLIS] skmm
template <typename T>
inline void skmm(side_t sidea,
                 uplo_t uploa,
                 conj_t conja,
                 trans_t transb,
                 dim_t m, dim_t n,
                 T alpha,
                 T *a, inc_t lda,
                 T *b, inc_t ldb,
                 T beta,
                 T *c, inc_t ldc);
template <> inline void skmm<float>(side_t sidea, uplo_t uploa, conj_t conja, trans_t transb, dim_t m, dim_t n, float alpha, float *a, inc_t lda, float *b, inc_t ldb, float beta, float *c, inc_t ldc)
{ bli_sskmm(sidea, uploa, conja, transb, m, n, &alpha, a, 1, lda, b, 1, ldb, &beta, c, 1, ldc); }
template <> inline void skmm<double>(side_t sidea, uplo_t uploa, conj_t conja, trans_t transb, dim_t m, dim_t n, double alpha, double *a, inc_t lda, double *b, inc_t ldb, double beta, double *c, inc_t ldc)
{ bli_dskmm(sidea, uploa, conja, transb, m, n, &alpha, a, 1, lda, b, 1, ldb, &beta, c, 1, ldc); }
template <> inline void skmm<ccscmplx>(side_t sidea, uplo_t uploa, conj_t conja, trans_t transb, dim_t m, dim_t n, ccscmplx alpha, ccscmplx *a, inc_t lda, ccscmplx *b, inc_t ldb, ccscmplx beta, ccscmplx *c, inc_t ldc)
{ bli_cskmm(sidea, uploa, conja, transb, m, n, (scomplex *)&alpha, (scomplex *)a, 1, lda, (scomplex *)b, 1, ldb, (scomplex *)&beta, (scomplex *)c, 1, ldc); }
template <> inline void skmm<ccdcmplx>(side_t sidea, uplo_t uploa, conj_t conja, trans_t transb, dim_t m, dim_t n, ccdcmplx alpha, ccdcmplx *a, inc_t lda, ccdcmplx *b, inc_t ldb, ccdcmplx beta, ccdcmplx *c, inc_t ldc)
{ bli_zskmm(sidea, uploa, conja, transb, m, n, (dcomplex *)&alpha, (dcomplex *)a, 1, lda, (dcomplex *)b, 1, ldb, (dcomplex *)&beta, (dcomplex *)c, 1, ldc); }

// [BLIS] skr2k
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

