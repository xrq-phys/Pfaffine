/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#include "blis.h"
#include "blalink.hh"
#include "colmaj.tcc"
#include "sktdi_diag.tcc"
#include <iostream>

template <typename T>
signed sktdi(uplo_t uplo,
             dim_t n,
             T *_A, inc_t ldA,
             T *_G, inc_t ldG,
             signed *iPov,
             T *_Work, dim_t lWork)
{
    // Define matrices.
    colmaj<T> A(_A, ldA);
    colmaj<T> G(_G, ldG);
    // Only one-line buffer is required.
    T *vA = _Work;

    // Inverse the central tri-diagonal part.
    sktdi_diag<T>(n, &A(0, 0), ldA);

    switch (uplo) {
    case BLIS_UPPER:
        for (int istep = n-3; istep >= 0; --istep) {
            // vA = invT * vG;
            // gemv<T>(BLIS_NO_TRANSPOSE, n, n, 1.0,
            //         &A(0, 0), ldA, &G(0, istep), 1, 0.0, vA, 1);
            gemv<T>(BLIS_NO_TRANSPOSE, n, n-istep-2, 1.0,
                    &A(0, istep+2), ldA, &G(istep+2, istep), 1, 0.0, vA, 1);

            // vA[:, istep+1] -= A * vG
            // vA[istep+1, :] += A * vG
            for (dim_t i = 0; i < n; ++i)
                A(i, istep+1) -= vA[i];
            for (dim_t i = 0; i < n; ++i)
                A(istep+1, i) += vA[i];
        }
        break;

    default:
        std::cerr << "SKTDI: Lower triangular storage not implemented. Sorry." << std::endl;
        return err_info(Pfaffine_NOT_IMPLEMNTED, 0);
        break;
    }

    // Check permutation and swap back.
    for (dim_t s = 0; s < n; ++s)
        if (iPov[s] != s) {
            dim_t t = s;
            for ( ; t < n; ++t)
                if (iPov[t] == s)
                    break;
            if (iPov[t] != s) {
                std::cerr << "SKTDI: iPov does not represent a valid permutation." << std::endl;
                return err_info(Pfaffine_BAD_REPRESENTATION, 7);
            }

            // Swap the data.
            swap<T>(n, &A(0, s), 1, &A(0, t), 1);
            swap<T>(n, &A(s, 0), ldA, &A(t, 0), ldA);
            // Swap the index.
            iPov[t] = iPov[s];
            iPov[s] = s;
        }

    return 0;
}
