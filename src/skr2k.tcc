/**
 * \file skr2k.tcc
 * Dispatcher and minimal implementation for SKR2k.
 * Minimal implementation is used when M is extremely small.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#include "blalink.hh"
#include "colmaj.tcc"
#include <iostream>

// Smaller than this size then minimal kernel will be called.
const dim_t ukr_dim_limit = 48;

template <typename T>
inline void ccukr_skr2k(uplo_t uploc,
                        trans_t transab,
                        dim_t m, dim_t k,
                        T alpha,
                        T *_A, inc_t ldA,
                        T *_B, inc_t ldB,
                        T beta,
                        T *_C, inc_t ldC)
{
    using namespace std;
    colmaj<T> A(_A, ldA);
    colmaj<T> B(_B, ldB);
    colmaj<T> C(_C, ldC);
    T one(1.0);

    if (transab == BLIS_NO_TRANSPOSE) {
        if (uploc == BLIS_UPPER) {
            for (signed j = 0; j < m; ++j) {
                if (beta != one) // 1.0 is exact.
                    for (signed i = 0; i < j; ++i)
                        C(i, j) *= beta;

                for (signed l = 0; l < k; ++l) {
                    T Ajl = alpha * A(j, l);
                    T Bjl = alpha * B(j, l);
                    for (signed i = 0; i < j; ++i)
                        C(i, j) += A(i, l) * Bjl - B(i, l) * Ajl;
                }
            }
        } else {
            cerr << "error: This Ukr is not implemented. Doing nothing." << endl;
        }
    } else {
        cerr << "error: This Ukr is not implemented. Doing nothing." << endl;
    }
}

template <typename T>
inline void skr2k(uplo_t uploc,
                  trans_t transab,
                  dim_t m, dim_t k,
                  T alpha,
                  T *_A, inc_t ldA,
                  T *_B, inc_t ldB,
                  T beta,
                  T *_C, inc_t ldC)
{
    if (transab == BLIS_NO_TRANSPOSE && uploc == BLIS_UPPER && m < ukr_dim_limit)
        ccukr_skr2k(uploc, transab, m, k, alpha, _A, ldA, _B, ldB, beta, _C, ldC);
    else
        ccbli_skr2k(uploc, transab, m, k, alpha, _A, ldA, _B, ldB, beta, _C, ldC);
}
