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
#include <memory>
#include <cstring>
#include "kersel.hh"
#include "blalink.hh"

// Macros for first-index-runs-fastest.
#define    A(i,j)    A[ (i) + (j)*(ldA) ]
#define    B(i,j)    B[ (i) + (j)*(ldB) ]
#define    C(i,j)    C[ (i) + (j)*(ldC) ]
#define pakA(i,j) pakA[ (i) + (j)*(mr) ]
#define pakB(i,j) pakB[ (i) + (j)*(nr) ]

// Block size of k.
const unsigned tracblk = 8;

template<typename T>
void uskr2k(unsigned n, unsigned k, T alpha, T *A, unsigned ldA, T *B, unsigned ldB, T beta, T *C, unsigned ldC,
            unsigned mr, unsigned nr, T *buffer)
{
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

    // Panel sizes
    std::size_t spaceAsz = mr * tracblk * sizeof(T) + 32;
    std::size_t spaceBsz = nr * tracblk * sizeof(T) + 32;
    // Allocate panels.
    void *spaceA = (void *)buffer;
    void *spaceB = (void *)((unsigned long)buffer + spaceAsz);
    T *pakA = (T *)std::align(32, mr * tracblk * sizeof(T), spaceA, spaceAsz);
    T *pakB = (T *)std::align(32, nr * tracblk * sizeof(T), spaceB, spaceBsz);

    // Number 1.
    T one = T(1.0);
    // minus alpha.
    T malpha = -alpha;

    for (unsigned uj = 0; uj < nblk; ++uj) {
        unsigned lenj = (uj+1==nblk) ? nblk_ : nr;
        for (unsigned ul = 0; ul < kblk; ++ul) {
            unsigned lenk = (ul+1==kblk) ? kblk_ : tracblk;
            T beta_ = (ul==1) ? beta : one;
            // Pack B into panels.
            // if (uj + 1 != nblk)
            for (unsigned l = 0; l < lenk; ++l)
                // TODO: change to more sophisticated routines.
                memcpy(&pakB(0, l), &B(uj*nr, l + ul*tracblk), sizeof(T) * lenj);
            for (unsigned ui = 0; ui < mblk; ++ui) {
                unsigned leni = (ui+1==mblk) ? mblk_ : mr;
                // Pack A into panels.
                for (unsigned l = 0; l < lenk; ++l)
                    memcpy(&pakA(0, l), &A(ui*mr, l + ul*tracblk), sizeof(T) * leni);
                // Starting indices.
                unsigned ist = ui*mr;
                unsigned jst = uj*nr;
                // TODO: prefetching.
                // Pick microkernel.
                if (ist + leni < jst)
                    if (mker_available<T>() && leni == mr && lenj == nr)
                        ugemmn(lenk, &alpha, pakA, pakB, &beta_, &C(ist, jst), ldC);
                    else
                        // Vanilla microkernel at off-diagonal.
                        for (unsigned j = 0; j < lenj; ++j)
                            for (unsigned l = 0; l < lenk; ++l) {
                                T bjl = pakB(j, l) * alpha;
                                for (unsigned i = 0; i < leni; ++i)
                                    C(ist + i, jst + j) =
                                        C(ist + i, jst + j) * beta_ + pakA(i, l) * bjl;
                            }
                else if (ist > jst + lenj)
                    if (mker_available<T>() && leni == mr && lenj == nr)
                        ugemmt(lenk, &malpha, pakB, pakA, &one, &C(jst, ist), ldC);
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
                        ugemmn(lenk, &alpha, pakA, pakB, &beta_, &C(ist, jst), ldC);
                        ugemmt(lenk, &malpha, pakB, pakA, &one, &C(jst, ist), ldC);
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
void skr2k(char uplo, char trans, unsigned n, unsigned k,
           T alpha, T *A, unsigned ldA, T *B, unsigned ldB, T beta, T *C, unsigned ldC, T *buffer)
{
    // Size to call directly interface GEMM.
    const unsigned mblk = 96;
    // Size of microblocks.
    unsigned mr, nr;
    set_blk_size<T>(&mr, &nr);

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
    // TODO: more checks.
    // TODO: special case when alpha==0

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

    for (unsigned mj = 0; mj < nblk; ++mj) {
        unsigned lenj = (mj+1 == nblk) ? nblk_ : mblk;
        for (unsigned mi = 0; mi < nblk; ++mi) {
            unsigned leni = (mi+1 == nblk) ? nblk_ : mblk;
            if (mi == mj)
                // SKR2K kernel.
                uskr2k<T>(lenj, k, alpha,
                          &A(mi*mblk, 0), ldA,
                          &B(mj*mblk, 0), ldB, beta,
                          &C(mi*mblk, mj*mblk), ldC, mr, nr, buffer);
            else if (mi < mj)
                // GEMM big-kernel.
                gemm('N', 'T', leni, lenj, k, alpha,
                     &A(mi*mblk, 0), ldA,
                     &B(mj*mblk, 0), ldB, beta,
                     &C(mi*mblk, mj*mblk), ldC);
            else
                // GEMM negative big-kernel.
                gemm('N', 'T', lenj, leni, k, -alpha,
                     &B(mj*mblk, 0), ldB,
                     &A(mi*mblk, 0), ldA, beta,
                     &C(mj*mblk, mi*mblk), ldC);
        }
    }
}


