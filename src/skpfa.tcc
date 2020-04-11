/**
 * \file skpfa.cc
 * For a Parlett-Reid method for calculation Pfaffian and inverse.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#include <iostream>
#include <blis/blis.h>

// Macros for first-index-runs-fastest.
#define   A(i,j)   A[ (i) + (j)*(ldA) ]
#define Sp1(i,j) Sp1[ (i) + (j)*(n)   ]
#define Sp2(i,j) Sp2[ (i) + (j)*(n)   ]
#define Sp3(i,j) Sp3[ (i) + (j)*(n)   ]
#define Sp4(i,j) Sp4[ (i) + (j)*(n)   ]
#define Sp5(i,j) Sp5[ (i) + (j)*(n)   ]

/**
 * \brief Calculate Pfaffian and optionally inverse of antisymmetric matrix A.
 *
 * \param uplo Only 'U' allowed.
 * \param n   Dimension of A.
 * \param A   Array (reference i.e. address of) A.
 * \param ldA Leading dimension size of A (row-skip).
 * \param inv Whether to compute inverse.
 * \param Sp[1-2] Two n*(npanel+1) scratchpad arrays.
 * \param Sp[3-4] Two n*n scratchpads for computing inverse, used only if inv=true.
 * \param Sp5 n*(npanel+1) scratchpad for computing inverse, untouched if inv=0.
 */
template <typename T>
T skpfa(char uplo, unsigned n, 
        T *A, unsigned ldA, unsigned inv,
        T *Sp1, T *Sp2, T *Sp3, T *Sp4, T *Sp5, unsigned npanel)
{
    // Error exit for not implemented.
    if (uplo != 'U' && uplo != 'u') {
        std::cerr << "Lower triangular storage not implemented. Sorry." << std::endl;
        std::_Exit(EXIT_FAILURE);
    }
    // Initialize updating-vector's scratchpad.
    T *vA = Sp1 + n*npanel;
    T *vG = Sp2 + n*npanel; // alpha_k: Gaussian elimination vector.
    T *vM, *kM;
    if (inv) { // Initialize M as identity.
        vM = Sp5 + n*npanel;
        kM = Sp4; // Reusable space.
        for (unsigned j = 0; j < n; ++j) {
            for (unsigned i = 0; i < n; ++i)
                Sp3(i, j) = 0.0;
            Sp3(j, j) = 1.0;
        }
    }

#ifdef _PR_Simple
    // Simple variant (without blocking or pivoting) for debugging.
    // vA[:] = A[:, 0];
    skslc<T>(n, 0, vA, A, ldA);
    for (unsigned istep = 0; istep < n-2; ++istep) {
        // k = istep+1, operating on column istep.
        // A(:, 0), kept from previous step.
        memcpy(vG, vA, n*sizeof(T));

        // vA = A[:, istep+1]
        skslc<T>(n, istep+1, vA, A, ldA);

        // Divide by alpha_k
        T alpha_k = vG[istep+1];
        for (unsigned i = 0; i < n; ++i)
            vG[i] /= alpha_k;
        for (unsigned i = 0; i <= istep+1; ++i)
            vG[i] = 0.0;

        // Here I'm not implementing additional skr2.
        // Directly call special case of skr2k<T>.
        skr2k<T>(uplo, 'N', n, 1, 1.0, vG, n, vA, n, 1.0, A, ldA);

#ifdef _Pfaff_Debug
        printf("After %d changes A=\n", istep+1);
        for (unsigned i = 0; i < n; ++i) {
            for (unsigned j = 0; j < n; ++j)
                printf("%12f ", A(i, j));
            putchar('\n');
        }
#endif

        // Update for inverse.
        if (inv) {
            memcpy(vM, &Sp3(0, istep+1), n*sizeof(T));
            ger(BLIS_NO_CONJUGATE, BLIS_NO_CONJUGATE, n, n, (T)-1.0, vM, 1, vG, 1, Sp3, 1, n);
        }
    }

#else
    for (unsigned ist = 0; ist < n-2; ist+=npanel) {
        unsigned lpanel = (ist+npanel >= n-2) ? n-2 - ist : npanel;
        // vA = A[:, ist], A's contribution to u_k.
        skslc<T>(n, ist, vA, A, ldA);
        for (unsigned i = 0; i < lpanel; ++i) {
            unsigned icur = ist + i;
            // Updated A[:, icur];
            memcpy(vG, vA, n*sizeof(T));
            
            // Pivoting
            if (fabs(icur) < 1e-3) {
                unsigned s = icur+1;
                unsigned t;
                T Gmax;
                findmax<T>(n, &t, &Gmax, vG);

                // Do swapping
                if (s < t) {
                    if (i != 0) {
                        // For unmerged updates, swap rows (col maj.).
                        swap(lpanel, &Sp1(s, 0), n, &Sp1(t, 0), n);
                        swap(lpanel, &Sp2(s, 0), n, &Sp2(t, 0), n);
                    }
                    // Swap A.
                    swap(s, &A(0, s), 1, &A(0, t), 1);
                    A(s, t) *= -1.0;
                    if (t > s+1) {
                        // TODO: This swap-and-flip needs independent kernel for max spd:
                        swap(t-s-1, &A(s+1, t), 1, &A(s, s+1), ldA);
                        for (unsigned j = s+1; j < t; ++j) {
                            A(j, t) *= -1;
                            A(s, j) *= -1;
                        }
                    }
                    if (t+1 < n)
                        swap(n-t-1, &A(s, t+1), ldA, &A(t, t+1), ldA);
                    // Update vectors.
                    vG[t] = vG[s];
                    vG[s] = Gmax;

                    // Inverse
                    if (inv) {
                        swap(n, &Sp3(0, s), 1, &Sp3(0, t), 1);
                        if (i != 0)
                            swap(lpanel, &Sp5(s, 0), n, &Sp5(t, 0), n);
                    }
                }
            }

            // vA = A[:, icur+1]
            skslc<T>(n, icur+1, vA, A, ldA);

            // Divide by alpha_k
            T alpha_k = vG[icur+1];
            for (unsigned j = 0; j < n; ++j)
                vG[j] /= alpha_k;
            for (unsigned j = 0; j <= icur+1; ++j)
                vG[j] = 0.0;

            // vA from Original A to updated components.
            if (i != 0) {
                gemv(BLIS_NO_TRANSPOSE, BLIS_NO_CONJUGATE, n, i,-1.0, Sp1, 1, n, &Sp2(icur+1, 0), n, 1.0, vA, 1);
                gemv(BLIS_NO_TRANSPOSE, BLIS_NO_CONJUGATE, n, i, 1.0, Sp2, 1, n, &Sp1(icur+1, 0), n, 1.0, vA, 1);
            }

            // M change.
            if (inv) {
                for (unsigned j = 0; j < i; ++j)
                    kM[j] = Sp5(icur+1, j);
                if (i != 0)
                    ger(BLIS_NO_CONJUGATE, BLIS_NO_CONJUGATE, n, i, (T)-1.0, vG, 1, kM, 1, Sp5, 1, n);
                memcpy(&Sp5(0, i), vG, n*sizeof(T));
            }

            // Write to change buffer.
            // TODO: I can directly do all the copying/updating to correct position.
            memcpy(&Sp1(0, i), vA, n*sizeof(T));
            memcpy(&Sp2(0, i), vG, n*sizeof(T));
        }

        // Apply transformation.
        skr2k<T>(uplo, 'N', n, lpanel, 1.0, Sp2, n, Sp1, n, 1.0, A, ldA);

#ifdef _Pfaff_Debug
        printf("After %d changes A=\n", ist+lpanel);
        for (unsigned i = 0; i < n; ++i) {
            for (unsigned j = 0; j < n; ++j)
                printf("%12f ", A(i, j));
            putchar('\n');
        }
#endif

        // Apply on M.
        if (inv) {
            // Borrows Sp4.
            memcpy(Sp4, &Sp3(0, ist+1), n*lpanel*sizeof(T));
            gemm(BLIS_NO_TRANSPOSE, BLIS_TRANSPOSE, 
                 n, n, lpanel, -1.0, Sp4, 1, n, Sp5, 1, n, 1.0, Sp3, 1, n);
        }
    }

#endif
    // Pfaffian
    T PfA = 1.0;
    for (unsigned i = 0; i < n-1; i+=2)
        PfA *= A(i, i+1);
    if (inv) {
        sktdi<T>(n, A, ldA);

        // M (A M'), i.e. M' A M in Wimmer's paper.
        // TODO: Write another skmm?
        gemm(BLIS_NO_TRANSPOSE, BLIS_TRANSPOSE,    n, n, n, (T)1.0, A, 1, ldA, Sp3, 1, n, (T)0.0, Sp4, 1, n);
        gemm(BLIS_NO_TRANSPOSE, BLIS_NO_TRANSPOSE, n, n, n, (T)1.0, Sp3, 1, n, Sp4, 1, n, (T)0.0, A, 1, ldA);
    }

    return std::abs(PfA); // Return only PfA (+ gauge), inverse is stored in A.
}

