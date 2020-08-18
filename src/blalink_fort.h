#include "stdlib.h"
#include "blis.h"


extern "C"
{
// converter from BLIS enum to BLAS char.
char trans2char(trans_t t)
{
    switch (t)
    {
    case BLIS_NO_TRANSPOSE:
        return 'N';

    case BLIS_TRANSPOSE:
        return 'T';

    case BLIS_CONJ_NO_TRANSPOSE:
        return 'C';

    default:
        abort();
        break;
    }
}


// gemm
void sgemm_(char *transa, char *transb, dim_t *m, dim_t *n, dim_t *k, float *alpha,
            float *a, dim_t *lda, float *b, dim_t *ldb, float *beta, float *c, dim_t *ldc);
void dgemm_(char *transa, char *transb, dim_t *m, dim_t *n, dim_t *k, double *alpha,
            double *a, dim_t *lda, double *b, dim_t *ldb, double *beta, double *c, dim_t *ldc);
void cgemm_(char *transa, char *transb, dim_t *m, dim_t *n, dim_t *k, void *alpha,
            void *a, dim_t *lda, void *b, dim_t *ldb, void *beta, void *c, dim_t *ldc);
void zgemm_(char *transa, char *transb, dim_t *m, dim_t *n, dim_t *k, void *alpha,
            void *a, dim_t *lda, void *b, dim_t *ldb, void *beta, void *c, dim_t *ldc);


// ger
void sger_(dim_t *m, dim_t *n, float *alpha,
           float *x, dim_t *incx, float *y, dim_t *incy, float *a, dim_t *lda);
void dger_(dim_t *m, dim_t *n, double *alpha,
           double *x, dim_t *incx, double *y, dim_t *incy, double *a, dim_t *lda);
void cgeru_(dim_t *m, dim_t *n, void *alpha,
            void *x, dim_t *incx, void *y, dim_t *incy, void *a, dim_t *lda);
void zgeru_(dim_t *m, dim_t *n, void *alpha,
            void *x, dim_t *incx, void *y, dim_t *incy, void *a, dim_t *lda);


// gemv
void sgemv_(char *trans, dim_t *m, dim_t *n, float *alpha, float *a, dim_t *lda,
            float *x, dim_t *incx, float *beta, float *y, dim_t *incy);
void dgemv_(char *trans, dim_t *m, dim_t *n, double *alpha, double *a, dim_t *lda,
            double *x, dim_t *incx, double *beta, double *y, dim_t *incy);
void cgemv_(char *trans, dim_t *m, dim_t *n, void *alpha, void *a, dim_t *lda,
            void *x, dim_t *incx, void *beta, void *y, dim_t *incy);
void zgemv_(char *trans, dim_t *m, dim_t *n, void *alpha, void *a, dim_t *lda,
            void *x, dim_t *incx, void *beta, void *y, dim_t *incy);


// swap
void sswap_(dim_t *n, float *sx, dim_t *incx, float *sy, dim_t *incy);
void dswap_(dim_t *n, double *sx, dim_t *incx, double *sy, dim_t *incy);
void cswap_(dim_t *n, void *sx, dim_t *incx, void *sy, dim_t *incy);
void zswap_(dim_t *n, void *sx, dim_t *incx, void *sy, dim_t *incy);

// axpy
void saxpy_(dim_t *n, float *alpha, float *sx, dim_t *incx, float *sy, dim_t *incy);
void daxpy_(dim_t *n, double *alpha, double *sx, dim_t *incx, double *sy, dim_t *incy);
void caxpy_(dim_t *n, void *alpha, void *sx, dim_t *incx, void *sy, dim_t *incy);
void zaxpy_(dim_t *n, void *alpha, void *sx, dim_t *incx, void *sy, dim_t *incy);

// dot
float sdot_(dim_t *n, float *sx, dim_t *incx, float *sy, dim_t *incy);
double ddot_(dim_t *n, double *sx, dim_t *incx, double *sy, dim_t *incy);
scomplex cdotc_(dim_t *n, void *sx, dim_t *incx, void *sy, dim_t *incy);
dcomplex zdotc_(dim_t *n, void *sx, dim_t *incx, void *sy, dim_t *incy);

}

