/**
 * \file blalink.hh
 * C++ to C type wrapper for BLAS.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#pragma once
#include <complex>
typedef std::complex<float>  scomplex;
typedef std::complex<double> dcomplex;
typedef struct {
    float real;
    float imag;
} cscomplex;
typedef struct {
    double real;
    double imag;
} cdcomplex;


// gemm
extern "C" {
void sgemm_(char *transa, char *transb, unsigned *m, unsigned *n, unsigned *k, float *alpha, 
            float *a, unsigned *lda, float *b, unsigned *ldb, float *beta, float *c, unsigned *ldc);
void dgemm_(char *transa, char *transb, unsigned *m, unsigned *n, unsigned *k, double *alpha, 
            double *a, unsigned *lda, double *b, unsigned *ldb, double *beta, double *c, unsigned *ldc);
void cgemm_(char *transa, char *transb, unsigned *m, unsigned *n, unsigned *k, void *alpha, 
            void *a, unsigned *lda, void *b, unsigned *ldb, void *beta, void *c, unsigned *ldc);
void zgemm_(char *transa, char *transb, unsigned *m, unsigned *n, unsigned *k, void *alpha, 
            void *a, unsigned *lda, void *b, unsigned *ldb, void *beta, void *c, unsigned *ldc);
}
inline void gemm(char transa, char transb, unsigned m, unsigned n, unsigned k, float alpha,
                 float *a, unsigned lda, float *b, unsigned ldb, float beta, float *c, unsigned ldc)
{ sgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); }
inline void gemm(char transa, char transb, unsigned m, unsigned n, unsigned k, double alpha,
                 double *a, unsigned lda, double *b, unsigned ldb, double beta, double *c, unsigned ldc)
{ dgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); }
inline void gemm(char transa, char transb, unsigned m, unsigned n, unsigned k, scomplex alpha,
                 scomplex *a, unsigned lda, scomplex *b, unsigned ldb, scomplex beta, scomplex *c, unsigned ldc)
{ cgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); }
inline void gemm(char transa, char transb, unsigned m, unsigned n, unsigned k, dcomplex alpha,
                 dcomplex *a, unsigned lda, dcomplex *b, unsigned ldb, dcomplex beta, dcomplex *c, unsigned ldc)
{ zgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); }


// ger
extern "C" {
void sger_(unsigned *m, unsigned *n, float *alpha, 
           float *x, unsigned *incx, float *y, unsigned *incy, float *a, unsigned *lda);
void dger_(unsigned *m, unsigned *n, double *alpha, 
           double *x, unsigned *incx, double *y, unsigned *incy, double *a, unsigned *lda);
void cgeru_(unsigned *m, unsigned *n, void *alpha, 
            void *x, unsigned *incx, void *y, unsigned *incy, void *a, unsigned *lda);
void zgeru_(unsigned *m, unsigned *n, void *alpha, 
            void *x, unsigned *incx, void *y, unsigned *incy, void *a, unsigned *lda);
}
inline void ger(unsigned m, unsigned n, float alpha, 
                float *x, unsigned incx, float *y, unsigned incy, float *a, unsigned lda)
{ sger_(&m, &n, &alpha, x, &incx, y, &incy, a, &lda); }
inline void ger(unsigned m, unsigned n, double alpha, 
                double *x, unsigned incx, double *y, unsigned incy, double *a, unsigned lda)
{ dger_(&m, &n, &alpha, x, &incx, y, &incy, a, &lda); }
inline void ger(unsigned m, unsigned n, scomplex alpha, 
                scomplex *x, unsigned incx, scomplex *y, unsigned incy, scomplex *a, unsigned lda)
{ cgeru_(&m, &n, &alpha, x, &incx, y, &incy, a, &lda); }
inline void ger(unsigned m, unsigned n, dcomplex alpha, 
                dcomplex *x, unsigned incx, dcomplex *y, unsigned incy, dcomplex *a, unsigned lda)
{ zgeru_(&m, &n, &alpha, x, &incx, y, &incy, a, &lda); }


// gemv
extern "C" {
void sgemv_(char *trans, unsigned *m, unsigned *n, float *alpha, float *a, unsigned *lda, 
            float *x, unsigned *incx, float *beta, float *y, unsigned *incy);
void dgemv_(char *trans, unsigned *m, unsigned *n, double *alpha, double *a, unsigned *lda, 
            double *x, unsigned *incx, double *beta, double *y, unsigned *incy);
void cgemv_(char *trans, unsigned *m, unsigned *n, void *alpha, void *a, unsigned *lda, 
            void *x, unsigned *incx, void *beta, void *y, unsigned *incy);
void zgemv_(char *trans, unsigned *m, unsigned *n, void *alpha, void *a, unsigned *lda, 
            void *x, unsigned *incx, void *beta, void *y, unsigned *incy);
}
inline void gemv(char trans, unsigned m, unsigned n, float alpha, float *a, unsigned lda,
                 float *x, unsigned incx, float beta, float *y, unsigned incy)
{ sgemv_(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy); }
inline void gemv(char trans, unsigned m, unsigned n, double alpha, double *a, unsigned lda,
                 double *x, unsigned incx, double beta, double *y, unsigned incy)
{ dgemv_(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy); }
inline void gemv(char trans, unsigned m, unsigned n, scomplex alpha, scomplex *a, unsigned lda,
                 scomplex *x, unsigned incx, scomplex beta, scomplex *y, unsigned incy)
{ cgemv_(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy); }
inline void gemv(char trans, unsigned m, unsigned n, dcomplex alpha, dcomplex *a, unsigned lda,
                 dcomplex *x, unsigned incx, dcomplex beta, dcomplex *y, unsigned incy)
{ zgemv_(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy); }


// swap
extern "C" {
void sswap_(unsigned *n, float *sx, unsigned *incx, float *sy, unsigned *incy);
void dswap_(unsigned *n, double *sx, unsigned *incx, double *sy, unsigned *incy);
void cswap_(unsigned *n, void *sx, unsigned *incx, void *sy, unsigned *incy);
void zswap_(unsigned *n, void *sx, unsigned *incx, void *sy, unsigned *incy);
}
inline void swap(unsigned n, float *x, unsigned incx, float *y, unsigned incy)
{ sswap_(&n, x, &incx, y, &incy); }
inline void swap(unsigned n, double *x, unsigned incx, double *y, unsigned incy)
{ dswap_(&n, x, &incx, y, &incy); }
inline void swap(unsigned n, scomplex *x, unsigned incx, scomplex *y, unsigned incy)
{ cswap_(&n, x, &incx, y, &incy); }
inline void swap(unsigned n, dcomplex *x, unsigned incx, dcomplex *y, unsigned incy)
{ zswap_(&n, x, &incx, y, &incy); }

// axpy
extern "C" {
void saxpy_(unsigned *n, float *alpha, float *sx, unsigned *incx, float *sy, unsigned *incy);
void daxpy_(unsigned *n, double *alpha, double *sx, unsigned *incx, double *sy, unsigned *incy);
void caxpy_(unsigned *n, void *alpha, void *sx, unsigned *incx, void *sy, unsigned *incy);
void zaxpy_(unsigned *n, void *alpha, void *sx, unsigned *incx, void *sy, unsigned *incy);
}
inline void axpy(unsigned n, float alpha, float *x, unsigned incx, float *y, unsigned incy)
{ saxpy_(&n, &alpha, x, &incx, y, &incy); }
inline void axpy(unsigned n, double alpha, double *x, unsigned incx, double *y, unsigned incy)
{ daxpy_(&n, &alpha, x, &incx, y, &incy); }
inline void axpy(unsigned n, scomplex alpha, scomplex *x, unsigned incx, scomplex *y, unsigned incy)
{ caxpy_(&n, &alpha, x, &incx, y, &incy); }
inline void axpy(unsigned n, dcomplex alpha, dcomplex *x, unsigned incx, dcomplex *y, unsigned incy)
{ zaxpy_(&n, &alpha, x, &incx, y, &incy); }

// dot
extern "C" {
float sdot_(unsigned *n, float *sx, unsigned *incx, float *sy, unsigned *incy);
double ddot_(unsigned *n, double *sx, unsigned *incx, double *sy, unsigned *incy);
cscomplex cdotc_(unsigned *n, void *sx, unsigned *incx, void *sy, unsigned *incy);
cdcomplex zdotc_(unsigned *n, void *sx, unsigned *incx, void *sy, unsigned *incy);
}
inline float dot(unsigned n, float *sx, unsigned incx, float *sy, unsigned incy)
{ return sdot_(&n, sx, &incx, sy, &incy); }
inline double dot(unsigned n, double *sx, unsigned incx, double *sy, unsigned incy)
{ return ddot_(&n, sx, &incx, sy, &incy); }
inline scomplex dot(unsigned n, scomplex *sx, unsigned incx, scomplex *sy, unsigned incy)
{ cscomplex inner = cdotc_(&n, sx, &incx, sy, &incy);
  return scomplex(inner.real, inner.imag); }
inline dcomplex dot(unsigned n, dcomplex *sx, unsigned incx, dcomplex *sy, unsigned incy)
{ cdcomplex inner = zdotc_(&n, sx, &incx, sy, &incy);
  return dcomplex(inner.real, inner.imag); }
