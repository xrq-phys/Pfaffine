/**
 * \file ogemm.tcc
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#include "colmaj.tcc"

/**
 * \brief Small-size 1-lv blocked GEMM program in 'N' - 'T' shape.
 */
template<typename T>
void ogemm(unsigned m, unsigned n, unsigned k,
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
    unsigned mblk = m / mr;
    unsigned nblk = n / nr;
    unsigned mblk_= m % mr;
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
    std::size_t spacePak  = (mr * mblk + nr * nblk) * k * sizeof(T);
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
    opack<T>(m, k, pakAsz, &A(), ldA, pakAbase, mr);
    opack<T>(n, k, pakBsz, &B(), ldB, pakBbase, nr);

    for (unsigned uj = 0; uj < nblk; ++uj) {
        unsigned lenj = (uj+1==nblk) ? nblk_ : nr;
        for (unsigned ul = 0; ul < kblk; ++ul) {
            unsigned lenk = (ul+1==kblk) ? kblk_ : tracblk;
            T beta_ = (ul==1) ? beta : one;
            // Selected packed B panels.
            colmaj<T> pakB(pakBbase + pakBsz * uj       // <packed 1st indices.
                                    + pakBpn * ul, nr); // <packed 2nd indices.
            for (unsigned ui = 0; ui < mblk; ++ui) {
                unsigned leni = (ui+1==mblk) ? mblk_ : mr;
                // Selected packed A panels.
                colmaj<T> pakA(pakAbase + pakAsz * ui       // <packed 1st indices.
                                        + pakApn * ul, mr); // <packed 2nd indices.
                // Starting indices.
                unsigned ist = ui*mr;
                unsigned jst = uj*nr;
                // Prefetch panels for next iteration.
                T *pakAnext,
                  *pakBnext;
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
            }
        }
    }
}

