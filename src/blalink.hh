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
// BLAS definitions
#ifdef BLAS_EXTERNAL
#include "blalink_fort.h"
#endif

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
#ifndef BLAS_EXTERNAL
#define BLALINK_MAC(cctype, ctype, cchar) \
    template <> inline void gemm<cctype>(trans_t transa, trans_t transb, \
                                         dim_t m, dim_t n, dim_t k, \
                                         cctype alpha, \
                                         cctype *a, inc_t lda, \
                                         cctype *b, inc_t ldb, \
                                         cctype beta, \
                                         cctype *c, inc_t ldc) \
    { \
        bli_##cchar##gemm(transa, transb, \
                          m, n, k, \
                          (ctype *)&alpha, \
                          (ctype *)a, 1, lda, \
                          (ctype *)b, 1, ldb, \
                          (ctype *)&beta, \
                          (ctype *)c, 1, ldc); \
    }
#else
#define BLALINK_MAC(cctype, ctype, cchar) \
    template <> inline void gemm<cctype>(trans_t transa, trans_t transb, \
                                         dim_t m, dim_t n, dim_t k, \
                                         cctype alpha, \
                                         cctype *a, inc_t lda, \
                                         cctype *b, inc_t ldb, \
                                         cctype beta, \
                                         cctype *c, inc_t ldc) \
    { \
        char ta = trans2char(transa), \
             tb = trans2char(transb); \
        cchar##gemm_(&ta, &tb, &m, &n, &k, \
                     (ctype *)&alpha, \
                     (ctype *)a, &lda, \
                     (ctype *)b, &ldb, \
                     (ctype *)&beta, \
                     (ctype *)c, &ldc); \
    }
#endif
BLALINK_MAC( float,    float,    s )
BLALINK_MAC( double,   double,   d )
BLALINK_MAC( ccscmplx, scomplex, c )
BLALINK_MAC( ccdcmplx, dcomplex, z )
#undef BLALINK_MAC


// ger
template <typename T>
inline void ger(dim_t m, dim_t n,
                T alpha,
                T *x, inc_t incx,
                T *y, inc_t incy,
                T *a, inc_t lda);
#ifndef BLAS_EXTERNAL
#define BLALINK_MAC(cctype, ctype, cchar, cfunc) \
    template <> inline void ger<cctype>(dim_t m, dim_t n, \
                                        cctype alpha, \
                                        cctype *x, inc_t incx, \
                                        cctype *y, inc_t incy, \
                                        cctype *a, inc_t lda) \
    { \
        bli_##cchar##ger(BLIS_NO_CONJUGATE, \
                         BLIS_NO_CONJUGATE, \
                         m, n, \
                         (ctype *)&alpha, \
                         (ctype *)x, incx, \
                         (ctype *)y, incy, \
                         (ctype *)a, 1, lda); \
    }
#else
#define BLALINK_MAC(cctype, ctype, cchar, cfunc) \
    template <> inline void ger<cctype>(dim_t m, dim_t n, \
                                        cctype alpha, \
                                        cctype *x, inc_t incx, \
                                        cctype *y, inc_t incy, \
                                        cctype *a, inc_t lda) \
    { \
        cchar##cfunc##_(&m, &n, \
                        (ctype *)&alpha, \
                        (ctype *)x, &incx, \
                        (ctype *)y, &incy, \
                        (ctype *)a, &lda); \
    }
#endif
BLALINK_MAC( float,    float,    s, ger  )
BLALINK_MAC( double,   double,   d, ger  )
BLALINK_MAC( ccscmplx, scomplex, c, geru )
BLALINK_MAC( ccdcmplx, dcomplex, z, geru )
#undef BLALINK_MAC


// gemv
template <typename T>
inline void gemv(trans_t trans,
                 dim_t m, dim_t n,
                 T alpha,
                 T *a, inc_t lda,
                 T *x, inc_t incx,
                 T beta,
                 T *y, inc_t incy);
#ifndef BLAS_EXTERNAL
#define BLALINK_MAC(cctype, ctype, cchar) \
    template <> inline void gemv<cctype>(trans_t trans, \
                                         dim_t m, dim_t n, \
                                         cctype alpha, \
                                         cctype *a, inc_t lda, \
                                         cctype *x, inc_t incx, \
                                         cctype beta, \
                                         cctype *y, inc_t incy) \
    { \
        bli_##cchar##gemv(trans, \
                          BLIS_NO_CONJUGATE, \
                          m, n, \
                          (ctype *)&alpha, \
                          (ctype *)a, 1, lda, \
                          (ctype *)x, incx, \
                          (ctype *)&beta, \
                          (ctype *)y, incy); \
    }
#else
#define BLALINK_MAC(cctype, ctype, cchar) \
    template <> inline void gemv<cctype>(trans_t trans, \
                                         dim_t m, dim_t n, \
                                         cctype alpha, \
                                         cctype *a, inc_t lda, \
                                         cctype *x, inc_t incx, \
                                         cctype beta, \
                                         cctype *y, inc_t incy) \
    { \
        char t = trans2char(trans); \
        cchar##gemv_(&t, &m, &n, \
                     (ctype *)&alpha, \
                     (ctype *)a, &lda, \
                     (ctype *)x, &incx, \
                     (ctype *)&beta, \
                     (ctype *)y, &incy); \
    }
#endif
BLALINK_MAC( float,    float,    s )
BLALINK_MAC( double,   double,   d )
BLALINK_MAC( ccscmplx, scomplex, c )
BLALINK_MAC( ccdcmplx, dcomplex, z )
#undef BLALINK_MAC


// trmm
template <typename T>
inline void trmm(side_t sidea, 
                 uplo_t uploa, 
                 trans_t transa, 
                 dim_t m, dim_t n,
                 T alpha, 
                 T *a, inc_t lda, 
                 T *b, inc_t ldb);
#ifndef BLAS_EXTERNAL
#define BLALINK_MAC(cctype, ctype, cchar) \
    template <> inline void trmm<cctype>(side_t sidea, uplo_t uploa, trans_t transa, \
                                         dim_t m, dim_t n, \
                                         cctype alpha, \
                                         cctype *a, inc_t lda, \
                                         cctype *b, inc_t ldb) \
    { \
        bli_##cchar##trmm(sidea, uploa, transa, \
                          BLIS_NONUNIT_DIAG, m, n, \
                          (ctype *)&alpha, \
                          (ctype *)a, 1, lda, \
                          (ctype *)b, 1, ldb); \
    }
#else
#define BLALINK_MAC(cctype, ctype, cchar) \
    template <> inline void trmm<cctype>(side_t sidea, uplo_t uploa, trans_t transa, \
                                         dim_t m, dim_t n, \
                                         cctype alpha, \
                                         cctype *a, inc_t lda, \
                                         cctype *b, inc_t ldb) \
    { \
        char ul = uplo2char(uploa), \
             si = side2char(sidea), \
             tr = trans2char(transa), \
             dg = 'N'; \
        cchar##trmm_(&si, &ul, &tr, &dg, \
                     &m, &n, \
                     (ctype *)&alpha, \
                     (ctype *)&a, &lda, \
                     (ctype *)&b, &ldb); \
    }
#endif
BLALINK_MAC( float,    float,    s )
BLALINK_MAC( double,   double,   d )
BLALINK_MAC( ccscmplx, scomplex, c )
BLALINK_MAC( ccdcmplx, dcomplex, z )
#undef BLALINK_MAC


// trmv
template <typename T>
inline void trmv(uplo_t uploa,
                 trans_t transa,
                 dim_t m,
                 T alpha,
                 T *a, inc_t lda,
                 T *x, inc_t incx);
#define BLALINK_MAC(cctype, ctype, cchar) \
    template <> inline void trmv<cctype>(uplo_t uploa, trans_t transa, dim_t m, \
                                         cctype alpha, \
                                         cctype *a, inc_t lda, \
                                         cctype *x, inc_t incx) \
    { \
        bli_##cchar##trmv(uploa, transa, \
                          BLIS_NONUNIT_DIAG, m, \
                          (ctype *)&alpha, \
                          (ctype *)a, 1, lda, \
                          (ctype *)x, incx); \
    }
/*
 * TRMV has alpha in its templated interface, which is absent in BLAS.
 *
#define BLALINK_MAC(cctype, ctype, cchar) \
    template <> inline void trmv<cctype>(uplo_t uploa, trans_t transa, dim_t m, \
                                         cctype alpha, \
                                         cctype *a, inc_t lda, \
                                         cctype *x, inc_t incx) \
    { \
        char ul = uplo2char(uploa), \
             tr = trans2char(transa), \
             dg = 'N'; \
        cchar##trmv_(&ul, &tr, &dg, &m, \
                     (ctype *)a, &lda, \
                     (ctype *)x, *incx); \
    }
 */
BLALINK_MAC( float,    float,    s )
BLALINK_MAC( double,   double,   d )
BLALINK_MAC( ccscmplx, scomplex, c )
BLALINK_MAC( ccdcmplx, dcomplex, z )
#undef BLALINK_MAC


// swap
template <typename T>
inline void swap(dim_t n,
                 T *x, inc_t incx,
                 T *y, inc_t incy);
#ifndef BLAS_EXTERNAL
#define BLALINK_MAC(cctype, ctype, cchar) \
    template <> inline void swap<cctype>(dim_t n, \
                                         cctype *x, inc_t incx, \
                                         cctype *y, inc_t incy) \
    { \
        bli_##cchar##swapv(n, \
                           (ctype *)x, incx, \
                           (ctype *)y, incy); \
    }
#else
#define BLALINK_MAC(cctype, ctype, cchar) \
    template <> inline void swap<cctype>(dim_t n, \
                                         cctype *x, inc_t incx, \
                                         cctype *y, inc_t incy) \
    { \
        cchar##swap_(&n, \
                     (ctype *)x, &incx, \
                     (ctype *)y, &incy); \
    }
#endif
BLALINK_MAC( float,    float,    s )
BLALINK_MAC( double,   double,   d )
BLALINK_MAC( ccscmplx, scomplex, c )
BLALINK_MAC( ccdcmplx, dcomplex, z )
#undef BLALINK_MAC


// axpy
template <typename T>
inline void axpy(dim_t n,
                 T alpha,
                 T *x, inc_t incx,
                 T *y, inc_t incy);
#ifndef BLAS_EXTERNAL
#define BLALINK_MAC(cctype, ctype, cchar) \
    template <> inline void axpy<cctype>(dim_t n, \
                                         cctype alpha, \
                                         cctype *x, inc_t incx, \
                                         cctype *y, inc_t incy) \
    { \
        bli_##cchar##axpyv(BLIS_NO_CONJUGATE, n, \
                           (ctype *)&alpha, \
                           (ctype *)x, incx, \
                           (ctype *)y, incy); \
    }
#else
#define BLALINK_MAC(cctype, ctype, cchar) \
    template <> inline void axpy<cctype>(dim_t n, \
                                         cctype alpha, \
                                         cctype *x, inc_t incx, \
                                         cctype *y, inc_t incy) \
    { \
        cchar##axpy_(&n, \
                     (ctype *)&alpha, \
                     (ctype *)x, &incx, \
                     (ctype *)y, &incy); \
    }
#endif
BLALINK_MAC( float,    float,    s )
BLALINK_MAC( double,   double,   d )
BLALINK_MAC( ccscmplx, scomplex, c )
BLALINK_MAC( ccdcmplx, dcomplex, z )
#undef BLALINK_MAC

// dot
template <typename T>
inline T dot(dim_t n,
             T *sx, inc_t incx,
             T *sy, inc_t incy);
#ifndef BLAS_EXTERNAL
#define BLALINK_MAC(cctype, ctype, cchar, cfunc) \
    template <> inline cctype dot<cctype>(dim_t n, \
                                          cctype *sx, inc_t incx, \
                                          cctype *sy, inc_t incy) \
    { \
        cctype rho; \
        bli_##cchar##dotv(BLIS_NO_CONJUGATE, \
                          BLIS_NO_CONJUGATE, \
                          n, \
                          (ctype *)sx, incx, \
                          (ctype *)sy, incy, \
                          (ctype *)&rho); \
        return rho; \
    }
#else
#define BLALINK_MAC(cctype, ctype, cchar, cfunc) \
    template <> inline cctype dot<cctype>(dim_t n, \
                                          cctype *sx, inc_t incx, \
                                          cctype *sy, inc_t incy) \
    { \
        ctype rho = cchar##cfunc##_(&n, \
                                    (ctype *)sx, &incx, \
                                    (ctype *)sy, &incy); \
        return *((cctype *)&rho); \
    }
#endif
BLALINK_MAC( float,    float,    s, dot  )
BLALINK_MAC( double,   double,   d, dot  )
BLALINK_MAC( ccscmplx, scomplex, c, dotc )
BLALINK_MAC( ccdcmplx, dcomplex, z, dotc )
#undef BLALINK_MAC


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
#define BLALINK_MAC(cctype, ctype, cchar) \
    template <> inline void skmm<cctype>(side_t sidea, uplo_t uploa, conj_t conja, trans_t transb, \
                                         dim_t m, dim_t n, \
                                         cctype alpha, \
                                         cctype *a, inc_t lda, \
                                         cctype *b, inc_t ldb, \
                                         cctype beta, \
                                         cctype *c, inc_t ldc) \
    { \
        bli_##cchar##skmm(sidea, uploa, conja, transb, \
                          m, n, \
                          (ctype *)&alpha, \
                          (ctype *)a, 1, lda, \
                          (ctype *)b, 1, ldb, \
                          (ctype *)&beta, \
                          (ctype *)c, 1, ldc); \
    }
BLALINK_MAC( float,    float,    s )
BLALINK_MAC( double,   double,   d )
BLALINK_MAC( ccscmplx, scomplex, c )
BLALINK_MAC( ccdcmplx, dcomplex, z )
#undef BLALINK_MAC


// [BLIS] skr2k
template <typename T>
inline void ccbli_skr2k(uplo_t uploc,
                        trans_t transab,
                        dim_t m, dim_t k,
                        T alpha,
                        T *a, inc_t lda,
                        T *b, inc_t ldb,
                        T beta,
                        T *c, inc_t ldc);
#define BLALINK_MAC(cctype, ctype, cchar) \
    template <> inline void ccbli_skr2k<cctype>(uplo_t uploc, trans_t transab, \
                                                dim_t m, dim_t k, \
                                                cctype alpha, \
                                                cctype *a, inc_t lda, \
                                                cctype *b, inc_t ldb, \
                                                cctype beta, \
                                                cctype *c, inc_t ldc) \
    { \
        bli_##cchar##skr2k(uploc, transab, transab, \
                           m, k, \
                           (ctype *)&alpha, \
                           (ctype *)a, 1, lda, \
                           (ctype *)b, 1, ldb, \
                           (ctype *)&beta, \
                           (ctype *)c, 1, ldc); \
    }
BLALINK_MAC( float,    float,    s )
BLALINK_MAC( double,   double,   d )
BLALINK_MAC( ccscmplx, scomplex, c )
BLALINK_MAC( ccdcmplx, dcomplex, z )
#undef BLALINK_MAC

