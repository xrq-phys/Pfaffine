/**
 * \file sktdf.tcc
 * Factorizes an even-rank skew-symmetric matrix into a tri-diagonal
 *   component and a series of Gaussian transformations.
 * It will be called Parlett-Reid decomposition later on.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#include "blis.h"
#include "blalink.hh"
#include "colmaj.tcc"
#include "findmax.tcc"
// SkBLAS with BLIS assembly or plain C++.
#include "skr2k.tcc"
#include "skslc.tcc"
#include "skswp.tcc"
#include <iostream>

template <typename T>
signed sktdf(uplo_t uplo,
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
    dim_t cflp = 0;

    // Initialize iPov space.
    for (dim_t i = 0; i < n; ++i)
        iPov[i] = i;

    // Define matrices.
    colmaj<T> A(_A, ldA);
    colmaj<T> G(_G, ldG);
    colmaj<T> Sp(_Work, n);

    // Column pointers.
    T *vA = nullptr;
    T *vG = nullptr;

    switch (uplo) {
    case BLIS_UPPER:
        for (dim_t ist = 0; ist < n-2; ist+=npanel) {
            dim_t lpanel = (ist+npanel >= n-2) ? n-2 - ist : npanel;
            // Location for saving transformations.
            colmaj<T> exG(&G(0, ist), G.ld);
            for (dim_t i = 0; i < lpanel; ++i) {
                dim_t icur = ist + i;
                vG = &exG(0, i);
                if (i == 0)
                    // vA = A[:, ist].
                    // skslc<T>(uplo, n, ist, vG, A, ldA);
                    // Copy lower-trangular non-zero components.
                    for (dim_t j = icur+1; j < n; ++j)
                        vG[j] = -A(ist, j);
                else
                    // Updated A[:, icur] is copied from previous itr.
                    // memcpy(vG, vA, n*sizeof(T));
                    // In fact full copying is not necessary. Only for lower-trangular parts.
                    memcpy(vG+icur, vA+icur, (n-icur)*sizeof(T));
                // Additional zeroizing required after simplification above.
                for (dim_t j = 0; j < icur+1; ++j)
                    vG[j] = 0.0;
                vA = &Sp(0, i);

                // Pivoting
                if (std::abs(vG[icur+1]) < 1e-3) {
                    inc_t s = icur+1;
                    inc_t t;
                    T Gmax;
                    findmax<T>(n, &t, &Gmax, vG, s);

                    // Do swapping.
                    if (s < t) {
                        // Sign-flip corresponding to the swap.
                        cflp++;
                        // Record permutation change.
                        dim_t itmp = iPov[s];
                        iPov[s] = iPov[t];
                        iPov[t] = itmp;

                        if (i != 0) {
                            // For unmerged updates, swap rows (col maj.).
                            swap<T>(i, &Sp(s, 0), Sp.ld, &Sp(t, 0), Sp.ld);
                        }
                        if (icur != 0)
                            // Swap the whole recorded history.
                            swap<T>(icur, &G(s, 0), G.ld, &G(t, 0), G.ld);
                        // Swap A.
                        skswp(uplo, n, &A(0, 0), A.ld, s, t);
                        // Update vectors.
                        vG[t] = vG[s];
                        vG[s] = Gmax;
                    }
                }

                // vA = A[:, icur+1]
                // skslc<T>(uplo, n, icur+1, vA, A, ldA);
                // Only upper-triagular is icur, icur+1
                for (dim_t j = 0; j < icur; ++j)
                    vA[j] = 0.0;
                vA[icur] = A(icur, icur+1);
                vA[icur+1] = 0.0;
                for (dim_t j = icur+2; j < n; ++j)
                    vA[j] = -A(icur+1, j);

                // Divide by alpha_k
                T alpha_k = vG[icur+1];
                for (dim_t j = icur+1; j < n; ++j)
                    vG[j] /= alpha_k;
                // for (unsigned j = 0; j <= icur+1; ++j)
                //     vG[j] = 0.0;
                // Only lower triangular is copied, zeroizing also reduced.
                // NB: if vG is fully copied at step-1, should do also vG[icur-1]=0.0;
                vG[icur+1] = 0.0;

                // Fast-update of A[:, icur:icur+lpanel].
                // Directly write to corresponding subdiagonal position of A.
                A(icur, icur+1) = -alpha_k;

                // vA from Original A to updated components, skipping zeros.
                if (i != 0) {
                    gemv<T>(BLIS_NO_TRANSPOSE, n-icur, i,-1.0, &Sp(icur, 0), Sp.ld, &exG(icur+1, 0), exG.ld, 1.0, vA+icur, 1);
                    gemv<T>(BLIS_NO_TRANSPOSE, n-icur, i, 1.0, &exG(icur, 0), exG.ld, &Sp(icur+1, 0), Sp.ld, 1.0, vA+icur, 1);
                }

                // Write to change buffer (already done).
                // memcpy(&Sp (0, i), vA, n*sizeof(T));
                // memcpy(&exG(0, i), vG, n*sizeof(T));
            }

            // Apply transformation.
            // skr2k<T>(BLIS_UPPER, BLIS_NO_TRANSPOSE, n, lpanel, 1.0,
            //          &exG(0, 0), exG.ld, &Sp(0, 0), Sp.ld, 1.0, &A(0, 0), A.ld);
            // Skip already cancelled by previous steps.
            // skr2k<T>(BLIS_UPPER, BLIS_NO_TRANSPOSE, n-ist, lpanel, 1.0,
            //          &exG(ist, 0), exG.ld, &Sp(ist, 0), Sp.ld, 1.0, &A(ist, ist), A.ld);
            // Rows & columns associated to this step are sure to be tridiagonal.
            // Skip these columns as well (by applying Fast-update to subdiagonals above).
            dim_t inext = ist + lpanel;
            if (n - inext > 1)
                skr2k<T>(BLIS_UPPER, BLIS_NO_TRANSPOSE, n-inext, lpanel, 1.0,
                         &exG(inext, 0), exG.ld, &Sp(inext, 0), Sp.ld, 1.0, &A(inext, inext), A.ld);

#ifdef _Pfaff_Debug
            printf("After %d changes A=\n", ist+lpanel);
            for (dim_t i = 0; i < n; ++i) {
                for (dim_t j = 0; j < n; ++j)
                    printf("%12f ", std::real(i < j ? i == j ? 0.0 : A(i, j) : -A(j, i)));
                putchar('\n');
            }
#endif
        }
        break;

    default:
        std::cerr << "SKTDF: Lower triangular storage not implemented. Sorry." << std::endl;
        return err_info(Pfaffine_NOT_IMPLEMNTED, 0);
    }
    // Save total number of flips into end of iPov.
    iPov[n] = cflp;

    return 0;
}
