/**
 * \file skr2k.cc
 * For a skew-symmetric matrix C, computes rank-2k update C = C*beta + alpha*(A B' - B A')
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#include "skr2k.hh"
#include <iostream>
#include <cstring>
#include <memory>
#include "kersel.hh"
#include "kersel_ext.hh"
#include "blalink.hh"
#include "thread.h"
#include "colmaj.tcc"
#include "floorsqrt.tcc"

// For microkernels.
#include "opack.tcc"
#include "ogemm.tcc"

// TODO: Rename and move this function to oskr2k.tcc.
template<typename T>
void uskr2k(unsigned n, unsigned k,
            T alpha,
            T *_A, unsigned ldA,
            T *_B, unsigned ldB,
            T beta,
            T *_C, unsigned ldC,
            unsigned mr, unsigned nr,
            T *buffer)
{
    colmaj<T> A(_A, ldA);
    colmaj<T> B(_B, ldB);
    colmaj<T> C(_C, ldC);
    // Determine microblock size.
    // Sorry for the confusion, but mbkl here shares meaning with nblk, not representing block size.
    unsigned mblk = n / mr;
    unsigned nblk = n / nr;
    unsigned mblk_= n % mr;
    unsigned nblk_= n % nr;
    unsigned kblk = k / tracblk;
    unsigned kblk_= k % tracblk;
    if (mblk_ > 0)
        ++mblk;
    else
        mblk_ = mr;
    if (nblk_ > 0)
        ++nblk;
    else
        nblk_ = nr;
    if (kblk_ > 0)
        ++kblk;
    else
        kblk_ = tracblk;

    // Number 1.
    T one = T(1.0);
    // minus alpha.
    T malpha = -alpha;

    // Align scratchpad space.
    std::size_t spacePak  = (mr * mblk + nr * nblk) * k;
    std::size_t spaceBuff = (mr * mblk + nr * nblk) * k * sizeof(T) + align_blk;
    void *buffer_ = buffer;
    T *buffPak = (T *)std::align(align_blk, spacePak, buffer_, spaceBuff);

    // Allocate panels.
    T *pakBbase = buffPak;
    T *pakAbase = buffPak + nr * nblk * k;

    // Pack the whole A & B.
    unsigned pakAsz = mr * k;
    unsigned pakBsz = nr * k;
    unsigned pakApn = mr * tracblk;
    unsigned pakBpn = nr * tracblk;
    opack<T>(n, k, pakAsz, &A(), ldA, pakAbase, mr);
    opack<T>(n, k, pakBsz, &B(), ldB, pakBbase, nr);

    for (unsigned uj = 0; uj < nblk; ++uj) {
        unsigned lenj = (uj+1==nblk) ? nblk_ : nr;
        for (unsigned ul = 0; ul < kblk; ++ul) {
            unsigned lenk = (ul+1==kblk) ? kblk_ : tracblk;
            T beta_ = (ul==1) ? beta : one;
            // Select packed B panels.
            colmaj<T> pakB(pakBbase + pakBsz * uj       // <packed 1st indices.
                                    + pakBpn * ul, nr); // <packed 2nd indices.
            for (unsigned ui = 0; ui < mblk; ++ui) {
                unsigned leni = (ui+1==mblk) ? mblk_ : mr;
                // Select packed A panels.
                colmaj<T> pakA(pakAbase + pakAsz * ui       // <packed 1st indices.
                                        + pakApn * ul, mr); // <packed 2nd indices.
                // Starting indices.
                unsigned ist = ui*mr;
                unsigned jst = uj*nr;
                // Prefetch panels for next iteration.
                T *pakAnext = nullptr,
                  *pakBnext = nullptr;
                if (ui+1 != mblk) {
                    pakBnext = &pakB();
                    pakAnext = pakAbase + pakAsz *(ui+1)
                                        + pakApn * ul;
                } else if (ul+1 != kblk) {
                    pakBnext = pakBbase + pakBsz * uj
                                        + pakBpn *(ul+1);
                    pakAnext = pakAbase + pakApn *(ul+1);
                } else {
                    pakBnext = pakBbase + pakBsz *(uj+1);
                    pakAnext = pakAbase;
                }
                // Pick microkernel.
                if (ist + leni <= jst)
                    if (mker_available<T>() && leni == mr && lenj == nr)
                        ugemmn(lenk, &alpha, &pakA(), &pakB(), &beta_, &C(ist, jst), ldC, pakAnext, pakBnext);
                    else if (extker_available<T>())
                        ugemmext(leni, lenj, lenk, &alpha, &pakA(), mr, &pakB(), nr, &beta_, &C(ist, jst), ldC,
                                 pakAnext, pakBnext);
                    else
                        // Vanilla microkernel at off-diagonal.
                        for (unsigned j = 0; j < lenj; ++j)
                            for (unsigned l = 0; l < lenk; ++l) {
                                T bjl = pakB(j, l) * alpha;
                                for (unsigned i = 0; i < leni; ++i)
                                    C(ist + i, jst + j) =
                                        C(ist + i, jst + j) * beta_ + pakA(i, l) * bjl;
                            }
                else if (ist >= jst + lenj)
                    if (mker_available<T>() && leni == mr && lenj == nr)
                        ugemmt(lenk, &malpha, &pakB(), &pakA(), &one, &C(jst, ist), ldC, pakAnext, pakBnext);
                    else if (extker_available<T>())
                        ugemmext(lenj, leni, lenk, &malpha, &pakB(), nr, &pakA(), mr, &one, &C(jst, ist), ldC,
                                 pakAnext, pakBnext);
                    else
                        for (unsigned i = 0; i < leni; ++i)
                            for (unsigned l = 0; l < lenk; ++l) {
                                T ail = pakA(i, l) * alpha;
                                for (unsigned j = 0; j < lenj; ++j)
                                    C(jst + j, ist + i) =
                                        C(jst + j, ist + i) * one - pakB(j, l) * ail;
                            }
                else
                    // Special microkernel at the diagonal.
                    if (mker_available<T>() && leni == mr && lenj == nr) {
                        // Update both contributions at diagonal.
                        // NOTE: This will cause some of lower triangular part overriden.
                        //       Still, no effect on result.
                        ugemmn(lenk, &alpha, &pakA(), &pakB(), &beta_, &C(ist, jst), ldC, &pakB(), &pakA());
                        ugemmt(lenk, &malpha, &pakB(), &pakA(), &one, &C(jst, ist), ldC, pakAnext, pakBnext);
                    } else if (extker_available<T>()) {
                        ugemmext(leni, lenj, lenk, &alpha, &pakA(), mr, &pakB(), nr, &beta_, &C(ist, jst), ldC,
                                 &pakB(), &pakA());
                        ugemmext(lenj, leni, lenk, &malpha, &pakB(), nr, &pakA(), mr, &one, &C(jst, ist), ldC,
                                 pakAnext, pakBnext);
                    } else
                        if (leni == lenj && ist == jst)
                            // Symmetric diagonal case.
                            for (unsigned j = 0; j < lenj; ++j)
                                for (unsigned l = 0; l < lenk; ++l) {
                                    T bjl = pakB(j, l) * alpha;
                                    T ajl = pakA(j, l) * alpha;
                                    for (unsigned i = 0; i < j; ++i) {
                                        unsigned i_ = i + ist;
                                        unsigned j_ = j + jst;
                                        C(i_, j_) = C(i_, j_) * beta_ +
                                            pakA(i, l) * bjl - pakB(i, l) * ajl;
                                    }
                                }
                        else
                            for (unsigned j = 0; j < lenj; ++j)
                                for (unsigned l = 0; l < lenk; ++l) {
                                    T bjl = pakB(j, l) * alpha;
                                    for (unsigned i = 0; i < leni; ++i) {
                                        // Indices in the 'big' block.
                                        unsigned i_ = i + ist;
                                        unsigned j_ = j + jst;
                                        if (i_ < j_)
                                            C(i_, j_) = C(i_, j_) * beta_ + pakA(i, l) * bjl;
                                        else if (i_ > j_)
                                            C(j_, i_) = C(j_, i_) * one - pakA(i, l) * bjl;
                                        else
                                            C(i_, j_) = 0.0;
                                    }
                                }
            }
        }
    }
}

template<typename T>
void skr2k(char uplo, char trans,
           unsigned n, unsigned k,
           T alpha,
           T *_A, unsigned ldA,
           T *_B, unsigned ldB,
           T beta,
           T *_C, unsigned ldC,
           T *buffer)
{
    colmaj<T> A(_A, ldA);
    colmaj<T> B(_B, ldB);
    colmaj<T> C(_C, ldC);
    // Size to call directly interface GEMM.
    const unsigned mblk = extblk;
    // Size of microblocks.
    unsigned mr, nr;
    set_blk_size<T>(&mr, &nr);

    // Calculate again nbitspbla.
    // TODO: Avoid calculating again for neat code.
    size_t pakAsz = k * mr,
           pakBsz = k * nr;
    size_t mmicroblk = extblk / mr + ((extblk % mr) ? 1 : 0);
    size_t nmicroblk = extblk / nr + ((extblk % nr) ? 1 : 0);
    size_t nbitspbla = sizeof(T) * (nmicroblk * pakBsz + mmicroblk * pakAsz) + align_blk;

    // Lo is not implemented. Sorry.
    if (uplo != 'U' && uplo != 'u') {
        std::cerr << "Lo is not implemented. Sorry." << std::endl;
        std::_Exit(EXIT_FAILURE);
    }
    // Trans is not implemented. Sorry.
    if (trans != 'N' && trans != 'n') {
        std::cerr << "Trans is not implemented. Sorry." << std::endl;
        std::_Exit(EXIT_FAILURE);
    }
    if (std::real(alpha) == 0.0 && std::imag(alpha) == 0.0)
        return;
    // TODO: More checks.

    // Big-blocking scheme, only for n.
    unsigned nblk = n / mblk;
    unsigned nblk_= n % mblk;
    nblk  =!nblk_ ? nblk  : nblk + 1;
    nblk_ = nblk_ ? nblk_ : mblk;
    // Exception
    if (nblk == 0) {
        nblk = 1;
        nblk_= n;
    }

#pragma omp parallel default(shared)
    {
        // Scratchpad corresponding to each thread.
        unsigned mp_id = omp_get_thread_num(),
                 mp_stride = omp_get_num_threads();
        T *buff_thread = (T *)((uint8_t *)buffer + (size_t)(nbitspbla * mp_id));

        // Manual scheduling.
        for (unsigned mij = mp_id; mij < nblk * (nblk + 1) / 2; mij += mp_stride) {
            unsigned mi, mj;
            glotr2ij(mij, mi, mj);
            unsigned lenj = (mj+1 == nblk) ? nblk_ : mblk;
            unsigned leni = (mi+1 == nblk) ? nblk_ : mblk;
            if (mij + mp_stride < nblk * (nblk + 1) / 2) {
                // Prefetching
                unsigned mi_n, mj_n;
                glotr2ij(mij + mp_stride, mi_n, mj_n);
                // Prefetch A & B.
#if defined(_Neon) || defined(_SVE)
                void *nextAiSq = &A(mi_n*mblk, 0);
                void *nextBiSq = &B(mi_n*mblk, 0);
                void *nextAjSq = &A(mj_n*mblk, 0);
                void *nextBjSq = &B(mj_n*mblk, 0);
                __asm__ volatile (
                    " ldr  x0, %[nextAiSq] \r\n" \
                    " ldr  x1, %[nextBiSq] \r\n" \
                    " prfm PLDL2STRM, [x0] \r\n" \
                    " prfm PLDL2STRM, [x1] \r\n"
                    :
                    : [nextAiSq] "m" (nextAiSq), \
                      [nextBiSq] "m" (nextBiSq)
                    : "x0", "x1"
                );
                if (mj_n != mj)
                    __asm__ volatile (
                        " ldr  x0, %[nextAiSq] \r\n" \
                        " ldr  x1, %[nextBiSq] \r\n" \
                        " prfm PLDL2STRM, [x0] \r\n" \
                        " prfm PLDL2STRM, [x1] \r\n"
                        :
                        : [nextAiSq] "m" (nextAiSq), \
                          [nextBiSq] "m" (nextBiSq)
                        : "x0", "x1"
                    );
// #elif defined(_Sandy) || defined(_Haswell) || defined(_SkylakeX)
// TODO: Add this Intel-specific prefetching.
#elif defined(__GNUC__)
                __builtin_prefetch(&A(mi_n*mblk, 0));
                __builtin_prefetch(&B(mi_n*mblk, 0));
                if (mj_n != mj) {
                    __builtin_prefetch(&A(mj_n*mblk, 0));
                    __builtin_prefetch(&B(mj_n*mblk, 0));
                }
#else
                // ANSI C has no prefetching. Doing nothing.
#endif
            }
            if (mi == mj) {
                // SKR2K kernel.
                uskr2k<T>(lenj, k, alpha,
                          &A(mj*mblk, 0), ldA,
                          &B(mj*mblk, 0), ldB, beta,
                          &C(mj*mblk, mj*mblk), ldC, mr, nr, buff_thread);
            } else { // GEMM extended-kernel.
#ifndef _Use_OGemm
                gemm('N', 'T', leni, lenj, k, alpha,
                     &A(mi*mblk, 0), ldA,
                     &B(mj*mblk, 0), ldB, beta,
                     &C(mi*mblk, mj*mblk), ldC);
#else
                ogemm<T>(leni, lenj, k, alpha,
                         &A(mi*mblk, 0), ldA,
                         &B(mj*mblk, 0), ldB, beta,
                         &C(mi*mblk, mj*mblk), ldC, mr, nr, buff_thread);
#endif
                // GEMM negative big-kernel.
#ifndef _Use_OGemm
                gemm('N', 'T', leni, lenj, k, -alpha,
                     &B(mi*mblk, 0), ldB,
                     &A(mj*mblk, 0), ldA, beta,
                     &C(mi*mblk, mj*mblk), ldC);
#else
                ogemm<T>(leni, lenj, k, -alpha,
                         &B(mi*mblk, 0), ldB,
                         &A(mj*mblk, 0), ldA, beta,
                         &C(mi*mblk, mj*mblk), ldC, mr, nr, buff_thread);
#endif
            }
        }
    }
}


