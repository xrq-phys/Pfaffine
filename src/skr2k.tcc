/**
 * \file skr2k.cc
 * For a skew-symmetric matrix C, computes rank-2k update C = C*beta + alpha*(A B' - B A')
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#include "typeswitch.hh"
#include <iostream>

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
            unsigned mr, unsigned nr, void *ugemm, auxinfo_t *aux, cntx_t *cntx, T *buffer)
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

    // Allocate panels.
    T *pakA = buffer;
    T *pakB = buffer + mr * k;
    T *tmpC = new T[mr * k];
    // Number 1.
    T one = T(1.0);

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
                // Very rough prefetching scheme.
                // bli_auxinfo_set_next_b(pakB, aux);
                // bli_auxinfo_set_next_a(&A((ui+1)*mblk, ul*kblk), aux);
                // Pick microkernel.
                if (ist + leni < jst)
                    if (false) // (ui + 1 != mblk && uj + 1 != nblk)
                        call_ugemm(ugemm, lenk, alpha, pakA, pakB, beta_, &C(ist, jst), 1, ldC, aux, cntx);
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
                    if (false) // (ui + 1 != mblk && uj + 1 != nblk)
                        call_ugemm(ugemm, lenk, -alpha, pakA, pakB, one, &C(jst, ist), ldC, 1, aux, cntx);
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
    delete[] tmpC;
}

template<typename T>
void skr2k(char uplo, char trans, unsigned n, unsigned k,
           T alpha, T *A, unsigned ldA, T *B, unsigned ldB, T beta, T *C, unsigned ldC)
{
    // Set kernel information (except 'next' entries).
    // TODO: query this information at higher levels could be better.
    dim_t mr, nr;
    auxinfo_t aux;
    cntx_t *cntx = bli_gks_query_cntx();
    void *ugemm = get_l3uker(C, BLIS_GEMM_UKR, cntx);
    set_blk_size(C, &mr, &nr, cntx);
    set_blis_is(C, &aux);

    // Scratchpad space.
    // TODO: move it to a higher level.
    // TODO: support dynamic. Though no >8x8 is present at the moment
    //       (There will be. After SVE there WILL be.)
    T buffer[128];
    // Size to call directly interface GEMM.
    const unsigned mblk = 32;

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
                          &C(mi*mblk, mj*mblk), ldC, mr, nr, ugemm, &aux, cntx, buffer);
            else if (mi < mj)
                // GEMM big-kernel.
                gemm(BLIS_NO_TRANSPOSE, BLIS_TRANSPOSE, leni, lenj, k, alpha,
                     &A(mi*mblk, 0), 1, ldA,
                     &B(mj*mblk, 0), 1, ldB, beta,
                     &C(mi*mblk, mj*mblk), 1, ldC);
            else
                // GEMM minus big-kernel.
                gemm(BLIS_NO_TRANSPOSE, BLIS_TRANSPOSE, lenj, leni, k, -alpha,
                     &B(mj*mblk, 0), 1, ldB,
                     &A(mi*mblk, 0), 1, ldA, beta,
                     &C(mj*mblk, mi*mblk), 1, ldC);
        }
    }
}


