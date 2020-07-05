/**
 * \file skpfa.tcc
 * For a Parlett-Reid method for calculation Pfaffian and inverse.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#include "skpfa.hh"
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <complex>
#include "blalink.hh"
#include "skr2k.hh"
#include "kersel.hh"
#include "thread.h"

#include "gcd.tcc"
#include "findmax.tcc"
#include "skslc.tcc"
#include "sktdi.tcc"

// Macros for first-index-runs-fastest.
// TODO: switch to colmaj.tcc.
#define   A(i,j)   A[ (i) + (j)*(ldA) ]
#define Sp1(i,j) Sp1[ (i) + (j)*(n)   ]
#define Sp2(i,j) Sp2[ (i) + (j)*(n)   ]
#define Sp3(i,j) Sp3[ (i) + (j)*(n)   ]
#define Sp4(i,j) Sp4[ (i) + (j)*(n)   ]
#define Sp5(i,j) Sp5[ (i) + (j)*(n)   ]
#define exG(i,j) exG[ (i) + (j)*(n)   ]

/**
 * \brief Calculate Pfaffian and optionally inverse of antisymmetric matrix A.
 *
 * \param uplo Only 'U' allowed.
 * \param n   Dimension of A.
 * \param A   Array (reference i.e. address of) A.
 * \param ldA Leading dimension size of A (row-skip).
 * \param inv Whether to compute inverse.
 * \param Sp1 An n*npanel scratchpad array.
 * \param Sp2 An n*npanel scratchpad array, used if inv=false (Sp3 would be used otherwise).
 * \param Sp3 An n*n scratchpad for computing inverse, used only if inv=true.
 * \param Sp4 v0.2 compatibility. Not used.
 * \param Sp5 v0.2 compatibility. Not used.
 */
template <typename T>
T skpfa(char uplo, unsigned n,
        T *A, unsigned ldA, unsigned inv,
        T *Sp1, T *Sp2, T *Sp3, T *Sp4, T *Sp5, unsigned npanel)
{
    // Allocates packing space.
    // TODO: Try avoiding in-place allocation.
    // TODO: Or allowing a global in-place allocation.
    unsigned mr, nr;
    set_blk_size<T>(&mr, &nr);
#ifndef _Manual_SqBlk
    // Automatically determines a square block size.
    if (n >= 2 * lcm(mr, nr) ||
        n % lcm(mr, nr) < 16 ||
        n - n / lcm(mr, nr) * lcm(mr, nr) < 16) {
        unsigned sqblk_auto = lcm(mr, nr);
        // The external block should not be too small.
        while (sqblk_auto < 64 && sqblk_auto * 2 <= n)
            sqblk_auto *= 2;
        set_sqblk_size(sqblk_auto);
    } else
#if defined(_SkylakeX) || defined(_SVE)
        set_sqblk_size(128);
#else
        set_sqblk_size(64);
#endif
#endif
    unsigned pakAsz = npanel * mr;
    unsigned pakBsz = npanel * nr;
    // Pivoting ordering.
    // TODO: Allocate iPov externally.
    unsigned *iPov = nullptr;
    // A-packing needs full size for maximum performance.
    T *SpBla = nullptr;
    // For sign of Pfaffian.
    unsigned cflp = 0;

    // Allocate iPov.
    iPov = (unsigned *)malloc(sizeof(unsigned) * n);
    // Set initial permutation data.
    for (unsigned i = 0; i < n; ++i)
        iPov[i] = i;
    // Use malloc() as VLAs are somehow bad as the allocation is hard to check.
    // Final align_block * 2 is for alignment padding.
    size_t mmicroblk = extblk / mr + ((extblk % mr) ? 1 : 0);
    size_t nmicroblk = extblk / nr + ((extblk % nr) ? 1 : 0);
    size_t nbitspbla = sizeof(T) * (nmicroblk * pakBsz + mmicroblk * pakAsz) + align_blk;
    SpBla = (T *)malloc(nbitspbla * omp_get_max_threads());
    if (SpBla == nullptr) {
        std::cerr << "Unable to allocate memory-packing scratchpads." << std::endl;
        std::_Exit(EXIT_FAILURE);
    }

    // Error exit for not implemented.
    if (uplo != 'U' && uplo != 'u') {
        std::cerr << "Lower triangular storage not implemented. Sorry." << std::endl;
        std::_Exit(EXIT_FAILURE);
    }
    // Updating-vector's scratchpad.
    // Default value (1-line buffer) used only when _PR_Simple is enabled.
    T *vA = Sp1;
    T *vG = Sp2; // alpha_k: Gaussian elimination vector.
    T *vM = nullptr, *kM = nullptr;

#ifdef _PR_Simple
    // Simple variant (without blocking or pivoting) for debugging.
    // vA[:] = A[:, 0];
    skslc<T>(n, 0, vA, A, ldA);
    for (unsigned istep = 0; istep < n-2; ++istep) {
        // Computing inverse requires logging vG vectors.
        if (inv)
            vG = &Sp3(0, istep);

        // k = istep+1, operating on column istep.
        // A(:, 0), kept from previous step.
        memcpy(vG, vA, n*sizeof(T));

        if (true) {
            unsigned s = istep+1;
            int t;
            T Gmax;
            // Pivoting.
            findmax<T>(n, &t, &Gmax, vG, s);

            if (int(s) < t) {
                cflp++;
                // Record permutation.
                unsigned itmp = iPov[s];
                iPov[s] = iPov[t];
                iPov[t] = itmp;

                // Unmerged G.
                if (inv)
                    swap(istep, &Sp3(s, 0), n, &Sp3(t, 0), n);

                // Swap A.
                swap(s, &A(0, s), 1, &A(0, t), 1);
                A(s, t) *= -1.0;
                if (t > s+1) {
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
            }
        }

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
        skr2k<T>(uplo, 'N', n, 1, 1.0, vG, n, vA, n, 1.0, A, ldA, SpBla);

#ifdef _Pfaff_Debug
        printf("After %d changes A=\n", istep+1);
        for (unsigned i = 0; i < n; ++i) {
            for (unsigned j = 0; j < n; ++j)
                printf("%12f ", std::real(i < j ? i == j ? 0.0 : A(i, j) : -A(j, i)));
            putchar('\n');
        }
#endif
    }

#else
    for (unsigned ist = 0; ist < n-2; ist+=npanel) {
        unsigned lpanel = (ist+npanel >= n-2) ? n-2 - ist : npanel;
        T *exG = nullptr;
        if (!inv)
            exG = Sp2;
        else
            // Inversion requires saving transformations.
            exG = &Sp3(0, ist);

        for (unsigned i = 0; i < lpanel; ++i) {
            unsigned icur = ist + i;
            vG = &exG(0, i);
            if (i == 0)
                // vA = A[:, ist].
                // skslc<T>(n, ist, vG, A, ldA);
                // Copy lower-trangular non-zero components.
                for (unsigned j = icur+1; j < n; ++j)
                    vG[j] = -A(ist, j);
            else
                // Updated A[:, icur] is copied from previous itr.
                // memcpy(vG, vA, n*sizeof(T));
                // Infact full copying is not necessary. Only for lower-trangular parts.
                memcpy(vG+icur, vA+icur, (n-icur)*sizeof(T));
            // Additional zeroizing required after simplification above.
            for (unsigned j = 0; j < icur+1; ++j)
                vG[j] = 0.0;
            vA = &Sp1(0, i);

            // Pivoting
            if (std::abs(vG[icur+1]) < 1e-3) {
                unsigned s = icur+1;
                int t;
                T Gmax;
                findmax<T>(n, &t, &Gmax, vG, s);

                // Do swapping.
                if (int(s) < t) {
                    // Sign-flip corresponding to the swap.
                    cflp++;
                    // Record permutation change.
                    unsigned itmp = iPov[s];
                    iPov[s] = iPov[t];
                    iPov[t] = itmp;

                    if (i != 0) {
                        // For unmerged updates, swap rows (col maj.).
                        swap(i, &Sp1(s, 0), n, &Sp1(t, 0), n);
                        if (!inv)
                            // Swap only unmerged vG.
                            swap(i, &exG(s, 0), n, &exG(t, 0), n);
                    }
                    if (inv)
                        if (icur != 0)
                            // Swap the whole recorded history.
                            swap(icur, &Sp3(s, 0), n, &Sp3(t, 0), n);
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
                }
            }

            // vA = A[:, icur+1]
            // skslc<T>(n, icur+1, vA, A, ldA);
            // Only upper-triagular is icur, icur+1
            for (unsigned j = 0; j < icur; ++j)
                vA[j] = 0.0;
            vA[icur] = A(icur, icur+1);
            vA[icur+1] = 0.0;
            for (unsigned j = icur+2; j < n; ++j)
                vA[j] = -A(icur+1, j);

            // Divide by alpha_k
            T alpha_k = vG[icur+1];
            for (unsigned j = icur+1 /*0*/; j < n; ++j)
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
                gemv('N', n-icur, i,-1.0, &Sp1(icur, 0), n, &exG(icur+1, 0), n, 1.0, vA+icur, 1);
                gemv('N', n-icur, i, 1.0, &exG(icur, 0), n, &Sp1(icur+1, 0), n, 1.0, vA+icur, 1);
            }

            // Write to change buffer (already done).
            // memcpy(&Sp1(0, i), vA, n*sizeof(T));
            // memcpy(&exG(0, i), vG, n*sizeof(T));
            // if (inv) ...
        }

        // Apply transformation.
        // skr2k<T>(uplo, 'N', n, lpanel, 1.0, exG, n, Sp1, n, 1.0, A, ldA, SpBla);
        // Skip already cancelled by previous steps.
        // skr2k<T>(uplo, 'N', n-ist, lpanel, 1.0,
        //          &exG(ist, 0), n, &Sp1(ist, 0), n, 1.0, &A(ist, ist), ldA, SpBla);
        // Rows & columns associated to this step are sure to be tridiagonal.
        // Skip these columns as well (by applying Fast-update to subdiagonals above).
        unsigned inext = ist + lpanel;
        if (n - inext > 1)
            skr2k<T>(uplo, 'N', n-inext, lpanel, 1.0,
                     &exG(inext, 0), n, &Sp1(inext, 0), n, 1.0, &A(inext, inext), ldA, SpBla);

#ifdef _Pfaff_Debug
        printf("After %d changes A=\n", ist+lpanel);
        for (unsigned i = 0; i < n; ++i) {
            for (unsigned j = 0; j < n; ++j)
                printf("%12f ", std::real(i < j ? i == j ? 0.0 : A(i, j) : -A(j, i)));
            putchar('\n');
        }
#endif
    }

#endif
    // Free memory-packing scratchpad.
    free(SpBla);

    // Pfaffian
    T PfA = 1.0;
    for (unsigned i = 0; i < n-1; i+=2)
        PfA *= A(i, i+1);
    if (inv) {
        sktdi<T>(n, A, ldA);

#ifdef _PI_Simple
        for (int istep = n-3; istep >= 0; --istep) {
            // vA = invT * vG;
            // gemv('N', n, n, 1.0,
            //      &A(0, 0), ldA, &Sp3(0, istep), 1, 0.0, vA, 1);
            gemv('N', n, n-istep-2, 1.0,
                 &A(0, istep+2), ldA, &Sp3(istep+2, istep), 1, 0.0, vA, 1);

            // vA[:, istep+1] -= A * vG
            // vA[istep+1, :] += A * vG
            for (unsigned i = 0; i < n; ++i)
                A(i, istep+1) -= vA[i];
            for (unsigned i = 0; i < n; ++i)
                A(istep+1, i) += vA[i];
        }
#else
        for (int ist = ((n-2)/npanel - !((n-2)%npanel))*npanel;
             ist >= 0; ist -= npanel) {
            unsigned lpanel = (ist+npanel > n-2) ? n-2-ist : npanel;

            // Compute inverse update.
            // gemm('N', 'N', n, lpanel, n, 1.0,
            //      &A(0, 0), ldA, &Sp3(0, ist), n, 0.0, Sp1, n);
            gemm('N', 'N', n, lpanel, n-ist-2, 1.0,
                 &A(0, ist+2), ldA, &Sp3(ist+2, ist), n, 0.0, Sp1, n);
            for (int i = lpanel-1; i >= 0; --i) {
                // TODO: Use ger instead of axpy.
                for (int j = 0; j < i; ++j)
                    axpy(n, -Sp3(ist+i+1, ist+j), &Sp1(0, i), 1, &Sp1(0, j), 1);

                // Inner-product part.
                if (i != 0)
                    for (int k = i; k < lpanel; ++k)
                        // Sp1(ist+k+1, i-1) += dot(n, &Sp3(0, ist+i-1), 1, &Sp1(0, k), 1);
                        Sp1(ist+k+1, i-1) += dot(n-ist-i-1, &Sp3(ist+i+1, ist+i-1), 1,
                                                            &Sp1(ist+i+1, k), 1);
            }

            // Write to inv(T).
            for (unsigned j = 0; j < lpanel; ++j)
                for (unsigned i = 0; i < n; ++i)
                    A(i, ist+1+j) -= Sp1(i, j);
            for (unsigned j = 0; j < lpanel; ++j)
                for (unsigned i = 0; i < n; ++i)
                    A(ist+1+j, i) += Sp1(i, j);
        }
#endif
        // Check permutation and swap back.
        for (unsigned s = 0; s < n; ++s)
            if (iPov[s] != s) {
                unsigned t = s;
                for ( ; t < n; ++t)
                    if (iPov[t] == s)
                        break;
                if (iPov[t] != s) {
                    std::cerr << "Permulation buffer broken." << std::endl;
                    std::_Exit(EXIT_FAILURE);
                }

                // Swap the data.
                swap(n, &A(0, s), 1, &A(0, t), 1);
                swap(n, &A(s, 0), ldA, &A(t, 0), ldA);
                // Swap the index.
                iPov[t] = iPov[s];
                iPov[s] = s;
            }
    }
    // Free Pivots.
    free(iPov);

    // Return only PfA, inverse is stored in A.
    return (cflp % 2) ? -PfA : PfA;
}

