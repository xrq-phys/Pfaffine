/**
 * \file skswp.tcc
 * Swap 2 row-columns of an antisymmetric matrix.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#pragma once
#include "colmaj.tcc"
#include "blalink.hh"
#include <iostream>

template <typename T> 
inline signed skswp(uplo_t uploA,
                    dim_t m,
                    T *A_, inc_t ldA,
                    dim_t s, dim_t t) 
{
    colmaj<T> A(A_, ldA);
    // Swap s and t if t<s.
    if (s > t) {
        dim_t r = s;
        s = t;
        t = r;
    }

    switch (uploA) {
    case BLIS_UPPER:
        swap<T>(s, &A(0, s), 1, &A(0, t), 1);
        A(s, t) *= -1.0;
        if (t > s+1) {
            // TODO: This swap-and-flip needs independent kernel for max spd:
            swap<T>(t-s-1, &A(s+1, t), 1, &A(s, s+1), A.ld);
            for (dim_t j = s+1; j < t; ++j) {
                A(j, t) *= -1;
                A(s, j) *= -1;
            }
        }
        if (t+1 < m)
            swap<T>(m-t-1, &A(s, t+1), A.ld, &A(t, t+1), A.ld);
        return 0;

    default:
        std::cerr << "SKSWP: Lower triangular storage not implemented. Sorry." << std::endl;
        return err_info(Pfaffine_NOT_IMPLEMNTED, 0);
    }
}

