/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#include "blis.h"
#include "blalink.hh"
#include "colmaj.tcc"
#include "findmax.tcc"
#include "skslc.tcc"
#include <iostream>

template <typename T>
signed sktdf(uplo_t uplo,
             dim_t n,
             T *_A, inc_t ldA,
             T *_G, inc_t ldG,
             signed *iPov,
             T *_Work, dim_t lWork)
{
    dim_t cflp = 0;

    // Initialize iPov space.
    for (dim_t i = 0; i < n; ++i)
        iPov[i] = i;

    // Define matrices.
    colmaj<T> A(_A, ldA);
    colmaj<T> G(_G, ldG);

    // Column pointers.
    T *vA = _Work;
    T *vG;

    switch (uplo) {
    case BLIS_UPPER:
        // Simple variant (without blocking or pivoting) for debugging.
        // vA[:] = A[:, 0];
        skslc<T>(n, 0, vA, A, ldA);
        for (dim_t istep = 0; istep < n-2; ++istep) {
            vG = &G(0, istep);

            // k = istep+1, operating on column istep.
            // A(:, 0), kept from previous step.
            memcpy(vG, vA, n*sizeof(T));

            if (true) {
                int s = istep+1;
                int t;
                T Gmax;
                // Pivoting.
                findmax<T>(n, &t, &Gmax, vG, s);

                if (s < t) {
                    cflp++;
                    // Record permutation.
                    dim_t itmp = iPov[s];
                    iPov[s] = iPov[t];
                    iPov[t] = itmp;

                    // Previous G.
                    swap<T>(istep, &G(s, 0), n, &G(t, 0), n);

                    // Swap A.
                    swap<T>(s, &A(0, s), 1, &A(0, t), 1);
                    A(s, t) *= -1.0;
                    if (t > s+1) {
                        swap<T>(t-s-1, &A(s+1, t), 1, &A(s, s+1), ldA);
                        for (dim_t j = s+1; j < t; ++j) {
                            A(j, t) *= -1;
                            A(s, j) *= -1;
                        }
                    }
                    if (t+1 < n)
                        swap<T>(n-t-1, &A(s, t+1), ldA, &A(t, t+1), ldA);

                    // Update vectors.
                    vG[t] = vG[s];
                    vG[s] = Gmax;
                }
            }

            // vA = A[:, istep+1]
            skslc<T>(n, istep+1, vA, A, ldA);

            // Divide by alpha_k
            T alpha_k = vG[istep+1];
            for (dim_t i = 0; i < n; ++i)
                vG[i] /= alpha_k;
            for (dim_t i = 0; i <= istep+1; ++i)
                vG[i] = 0.0;

            // Here I'm not implementing additional skr2.
            // Directly call special case of skr2k<T>.
            skr2k<T>(BLIS_UPPER, BLIS_NO_TRANSPOSE, n, 1, 1.0, vG, ldG, vA, n, 1.0, A, ldA);

#ifdef _Pfaff_Debug
            printf("After %d changes A=\n", istep+1);
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
