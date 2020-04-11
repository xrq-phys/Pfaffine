/**
 * \file typeswitch.tcc
 * Mostly a C++ to C type wrapper for BLIS microkernels.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#include <blis/blis.h>
#include <complex>

// set_blk_size {
inline void set_blk_size(const float *T, dim_t *mr, dim_t *nr, cntx_t *cntx)
{ *mr = bli_cntx_get_blksz_def_dt(BLIS_FLOAT, BLIS_MR, cntx);
  *nr = bli_cntx_get_blksz_def_dt(BLIS_FLOAT, BLIS_NR, cntx); }
inline void set_blk_size(const double *T, dim_t *mr, dim_t *nr, cntx_t *cntx)
{ *mr = bli_cntx_get_blksz_def_dt(BLIS_DOUBLE, BLIS_MR, cntx);
  *nr = bli_cntx_get_blksz_def_dt(BLIS_DOUBLE, BLIS_NR, cntx); }
inline void set_blk_size(const std::complex<float> *T, dim_t *mr, dim_t *nr, cntx_t *cntx)
{ *mr = bli_cntx_get_blksz_def_dt(BLIS_SCOMPLEX, BLIS_MR, cntx);
  *nr = bli_cntx_get_blksz_def_dt(BLIS_SCOMPLEX, BLIS_NR, cntx); }
inline void set_blk_size(const std::complex<double> *T, dim_t *mr, dim_t *nr, cntx_t *cntx)
{ *mr = bli_cntx_get_blksz_def_dt(BLIS_DCOMPLEX, BLIS_MR, cntx);
  *nr = bli_cntx_get_blksz_def_dt(BLIS_DCOMPLEX, BLIS_NR, cntx); }
// }

// get_l3uker {
inline void *get_l3uker(const float *T, l3ukr_t type, cntx_t *cntx)
{ return bli_cntx_get_l3_nat_ukr_dt(BLIS_FLOAT, type, cntx); }
inline void *get_l3uker(const double *T, l3ukr_t type, cntx_t *cntx)
{ return bli_cntx_get_l3_nat_ukr_dt(BLIS_DOUBLE, type, cntx); }
inline void *get_l3uker(const std::complex<float> *T, l3ukr_t type, cntx_t *cntx)
{ return bli_cntx_get_l3_nat_ukr_dt(BLIS_SCOMPLEX, type, cntx); }
inline void *get_l3uker(const std::complex<double> *T, l3ukr_t type, cntx_t *cntx)
{ return bli_cntx_get_l3_nat_ukr_dt(BLIS_DCOMPLEX, type, cntx); }
// }

extern "C" { // call_ugemm (plain pointer)
void call_sugemm(sgemm_ukr_ft ugemm, dim_t k, float *alpha, float *a, float *b,
                 float *beta, float *c, inc_t rsc, inc_t csc, auxinfo_t *aux, cntx_t *cntx)
{ ugemm(k, alpha, a, b, beta, c, rsc, csc, aux, cntx); }
void call_dugemm(dgemm_ukr_ft ugemm, dim_t k, double *alpha, double *a, double *b,
                 double *beta, double *c, inc_t rsc, inc_t csc, auxinfo_t *aux, cntx_t *cntx)
{ ugemm(k, alpha, a, b, beta, c, rsc, csc, aux, cntx); }
void call_cugemm(cgemm_ukr_ft ugemm, dim_t k, scomplex *alpha, scomplex *a, scomplex *b,
                 scomplex *beta, scomplex *c, inc_t rsc, inc_t csc, auxinfo_t *aux, cntx_t *cntx)
{ ugemm(k, alpha, a, b, beta, c, rsc, csc, aux, cntx); }
void call_zugemm(zgemm_ukr_ft ugemm, dim_t k, dcomplex *alpha, dcomplex *a, dcomplex *b,
                 dcomplex *beta, dcomplex *c, inc_t rsc, inc_t csc, auxinfo_t *aux, cntx_t *cntx)
{ ugemm(k, alpha, a, b, beta, c, rsc, csc, aux, cntx); }
} // {
inline void call_ugemm(void *ugemm, dim_t k, float alpha, float *a, float *b,
                       float beta, float *c, inc_t rsc, inc_t csc, auxinfo_t *aux, cntx_t *cntx)
{ call_sugemm((sgemm_ukr_ft)ugemm, k, &alpha, a, b, &beta, c, rsc, csc, aux, cntx); }
inline void call_ugemm(void *ugemm, dim_t k, double alpha, double *a, double *b,
                       double beta, double *c, inc_t rsc, inc_t csc, auxinfo_t *aux, cntx_t *cntx)
{ call_dugemm((dgemm_ukr_ft)ugemm, k, &alpha, a, b, &beta, c, rsc, csc, aux, cntx); }
inline void call_ugemm(void *ugemm, dim_t k, std::complex<float> alpha, std::complex<float> *a, std::complex<float> *b,
                       std::complex<float> beta, std::complex<float> *c, inc_t rsc, inc_t csc, auxinfo_t *aux, cntx_t *cntx)
{ call_cugemm((cgemm_ukr_ft)ugemm, k, (scomplex *)&alpha, (scomplex *)a, (scomplex *)b, (scomplex *)&beta, (scomplex *)c, rsc, csc, aux, cntx); }
inline void call_ugemm(void *ugemm, dim_t k, std::complex<double> alpha, std::complex<double> *a, std::complex<double> *b,
                       std::complex<double> beta, std::complex<double> *c, inc_t rsc, inc_t csc, auxinfo_t *aux, cntx_t *cntx)
{ call_zugemm((zgemm_ukr_ft)ugemm, k, (dcomplex *)&alpha, (dcomplex *)a, (dcomplex *)b, (dcomplex *)&beta, (dcomplex *)c, rsc, csc, aux, cntx); }
// }

// set_blis_is {
inline void set_blis_is(const float *T, auxinfo_t *aux)
{ bli_auxinfo_set_is_a(1, aux);
  bli_auxinfo_set_is_b(1, aux); }
inline void set_blis_is(const double *T, auxinfo_t *aux)
{ bli_auxinfo_set_is_a(1, aux);
  bli_auxinfo_set_is_b(1, aux); }
// TODO: check correctness of complex is=2.
inline void set_blis_is(const std::complex<float> *T, auxinfo_t *aux)
{ bli_auxinfo_set_is_a(2, aux);
  bli_auxinfo_set_is_b(2, aux); }
inline void set_blis_is(const std::complex<double> *T, auxinfo_t *aux)
{ bli_auxinfo_set_is_a(2, aux);
  bli_auxinfo_set_is_b(2, aux); }
// }

// gemm {
inline void gemm(trans_t transa, trans_t transb, dim_t m, dim_t n, dim_t k, float alpha,
                 float *a, inc_t rsa, inc_t csa, float *b, inc_t rsb, inc_t csb, float beta,
                 float *c, inc_t rsc, inc_t csc)
{ bli_sgemm(transa, transb, m, n, k, &alpha, a, rsa, csa, b, rsb, csb,
            &beta, c, rsc, csc); }
inline void gemm(trans_t transa, trans_t transb, dim_t m, dim_t n, dim_t k, double alpha,
                 double *a, inc_t rsa, inc_t csa, double *b, inc_t rsb, inc_t csb, double beta,
                 double *c, inc_t rsc, inc_t csc)
{ bli_dgemm(transa, transb, m, n, k, &alpha, a, rsa, csa, b, rsb, csb,
            &beta, c, rsc, csc); }
inline void gemm(trans_t transa, trans_t transb, dim_t m, dim_t n, dim_t k, std::complex<float> alpha,
                 std::complex<float> *a, inc_t rsa, inc_t csa, std::complex<float> *b, inc_t rsb, inc_t csb,
                 std::complex<float> beta,
                 std::complex<float> *c, inc_t rsc, inc_t csc)
{ bli_cgemm(transa, transb, m, n, k, (scomplex *)&alpha, (scomplex *)a, rsa, csa, (scomplex *)b, rsb, csb,
            (scomplex *)&beta, (scomplex *)c, rsc, csc); }
inline void gemm(trans_t transa, trans_t transb, dim_t m, dim_t n, dim_t k, std::complex<double> alpha,
                 std::complex<double> *a, inc_t rsa, inc_t csa, std::complex<double> *b, inc_t rsb, inc_t csb,
                 std::complex<double> beta,
                 std::complex<double> *c, inc_t rsc, inc_t csc)
{ bli_zgemm(transa, transb, m, n, k, (dcomplex *)&alpha, (dcomplex *)a, rsa, csa, (dcomplex *)b, rsb, csb,
            (dcomplex *)&beta, (dcomplex *)c, rsc, csc); }
// }

// ger {
inline void ger(conj_t conjx, conj_t conjy, dim_t m, dim_t n, float alpha, float *x, inc_t incx, float *y, inc_t incy, 
                float *a, inc_t rsa, inc_t csa)
{ bli_sger(conjx, conjy, m, n, &alpha, x, incx, y, incy, a, rsa, csa); }
inline void ger(conj_t conjx, conj_t conjy, dim_t m, dim_t n, double alpha, double *x, inc_t incx, double *y, inc_t incy, 
                double *a, inc_t rsa, inc_t csa)
{ bli_dger(conjx, conjy, m, n, &alpha, x, incx, y, incy, a, rsa, csa); }
inline void ger(conj_t conjx, conj_t conjy, dim_t m, dim_t n, std::complex<float> alpha, 
                std::complex<float> *x, inc_t incx, std::complex<float> *y, inc_t incy, 
                std::complex<float> *a, inc_t rsa, inc_t csa)
{ bli_cger(conjx, conjy, m, n, (scomplex *)&alpha, (scomplex *)x, incx, (scomplex *)y, incy, (scomplex *)a, rsa, csa); }
inline void ger(conj_t conjx, conj_t conjy, dim_t m, dim_t n, std::complex<double> alpha, 
                std::complex<double> *x, inc_t incx, std::complex<double> *y, inc_t incy, 
                std::complex<double> *a, inc_t rsa, inc_t csa)
{ bli_zger(conjx, conjy, m, n, (dcomplex *)&alpha, (dcomplex *)x, incx, (dcomplex *)y, incy, (dcomplex *)a, rsa, csa); }
// }

// gemv {
inline void gemv(trans_t transa, conj_t conjx, dim_t m, dim_t n, float alpha, float *a, inc_t rsa, inc_t csa,
                 float *x, inc_t incx, float beta, float *y, inc_t incy)
{ bli_sgemv(transa, conjx, m, n, &alpha, a, rsa, csa, x, incx, &beta, y, incy); }
inline void gemv(trans_t transa, conj_t conjx, dim_t m, dim_t n, double alpha, double *a, inc_t rsa, inc_t csa,
                 double *x, inc_t incx, double beta, double *y, inc_t incy)
{ bli_dgemv(transa, conjx, m, n, &alpha, a, rsa, csa, x, incx, &beta, y, incy); }
inline void gemv(trans_t transa, conj_t conjx, dim_t m, dim_t n, std::complex<float> alpha, 
                 std::complex<float> *a, inc_t rsa, inc_t csa,
                 std::complex<float> *x, inc_t incx, std::complex<float> beta, std::complex<float> *y, inc_t incy)
{ bli_cgemv(transa, conjx, m, n, (scomplex *)&alpha, (scomplex *)a, rsa, csa, 
            (scomplex *)x, incx, (scomplex *)&beta, (scomplex *)y, incy); }
inline void gemv(trans_t transa, conj_t conjx, dim_t m, dim_t n, std::complex<double> alpha, 
                 std::complex<double> *a, inc_t rsa, inc_t csa,
                 std::complex<double> *x, inc_t incx, std::complex<double> beta, std::complex<double> *y, inc_t incy)
{ bli_zgemv(transa, conjx, m, n, (dcomplex *)&alpha, (dcomplex *)a, rsa, csa, 
            (dcomplex *)x, incx, (dcomplex *)&beta, (dcomplex *)y, incy); }
// }

// swap {
inline void swap(dim_t n, float *x, inc_t incx, float *y, inc_t incy)
{ bli_sswapv(n, x, incx, y, incy); }
inline void swap(dim_t n, double *x, inc_t incx, double *y, inc_t incy)
{ bli_dswapv(n, x, incx, y, incy); }
inline void swap(dim_t n, std::complex<float> *x, inc_t incx, std::complex<float> *y, inc_t incy)
{ bli_cswapv(n, (scomplex *)x, incx, (scomplex *)y, incy); }
inline void swap(dim_t n, std::complex<double> *x, inc_t incx, std::complex<double> *y, inc_t incy)
{ bli_zswapv(n, (dcomplex *)x, incx, (dcomplex *)y, incy); }
// }

