/**
 * \file sktdi.tcc
 * Inverts an even-rank skew-symmetric matrix from a Parlett-Reid
 *   decomposition (done by sktdf<T>).
 *
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
    // TODO: Determine automatically npanel.
    dim_t npanel = lWork / n;
    if (npanel > n-2)
        npanel = n-2;

    // Define matrices.
    colmaj<T> A(_A, ldA);
    colmaj<T> G(_G, ldG);
    colmaj<T> Sp(_Work, n);

    // Inverse the central tri-diagonal part.
    sktdi_diag<T>(n, &A(0, 0), ldA);

    switch (uplo) {
    case BLIS_UPPER:
        for (signed ist = ((n-2)/npanel - !((n-2)%npanel))*npanel;
             ist >= 0; ist -= npanel) {
            dim_t lpanel = (ist+npanel > n-2) ? n-2-ist : npanel;

            // Compute inverse update.
            // gemm<T>(BLIS_NO_TRANSPOSE, BLIS_NO_TRANSPOSE, n, lpanel, n, 1.0,
            //         &A(0, 0), ldA, &G(0, ist), ldG, 0.0, &Sp(0, 0), n);
            gemm<T>(BLIS_NO_TRANSPOSE, BLIS_NO_TRANSPOSE, n, lpanel, n-ist-2, 1.0,
                    &A(0, ist+2), ldA, &G(ist+2, ist), ldG, 0.0, &Sp(0, 0), n);
            for (int i = lpanel-1; i >= 0; --i) {
                for (int j = 0; j < i; ++j)
                    axpy(n, -G(ist+i+1, ist+j), &Sp(0, i), 1, &Sp(0, j), 1);

                // Inner-product part.
                if (i != 0)
                    for (int k = i; k < lpanel; ++k)
                        // Sp(ist+k+1, i-1) += dot(n, &G(0, ist+i-1), 1, &Sp(0, k), 1);
                        Sp(ist+k+1, i-1) += dot(n-ist-i-1, &G(ist+i+1, ist+i-1), 1,
                                                           &Sp(ist+i+1, k), 1);
            }

            // Write to inv(T).
            for (dim_t j = 0; j < lpanel; ++j)
                for (dim_t i = 0; i < n; ++i)
                    A(i, ist+1+j) -= Sp(i, j);
            for (dim_t j = 0; j < lpanel; ++j)
                for (dim_t i = 0; i < n; ++i)
                    A(ist+1+j, i) += Sp(i, j);
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