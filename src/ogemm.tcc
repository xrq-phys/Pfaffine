/**
 * \file ogemm.tcc
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

/**
 * \brief Small-size 1-lv blocked GEMM program in 'N' - 'T' shape.
 */
template<typename T>
void ogemm(unsigned m, unsigned n, unsigned k, 
           T alpha, T *A, unsigned ldA, T *B, unsigned ldB, 
           T beta, T *C, unsigned ldC, unsigned mr, unsigned nr, T *buffer)
{
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

    // Panel sizes.
    std::size_t spaceBsz = align_blk + nr * tracblk  * sizeof(T);
    std::size_t spaceAsz = align_blk + mr * mblk * k * sizeof(T);

    // Allocate panels.
    void *spaceB = (void *)buffer;
    void *spaceA = (void *)((uint8_t *)buffer + spaceBsz);
    T *pakB     = (T *)std::align(align_blk, nr * tracblk,  spaceB, spaceBsz);
    T *pakAbase = (T *)std::align(align_blk, mr * mblk * k, spaceA, spaceAsz);

    // Pack the whole A.
    unsigned pakAsz = mr * k;
    unsigned pakApn = mr * tracblk;
    for (unsigned ui = 0; ui < mblk; ++ui) {
        unsigned leni = (ui+1==mblk) ? mblk_ : mr;
        T *pakA = pakAbase + pakAsz * ui;
        for (unsigned l = 0; l < k; ++l)
            memcpy(&pakA(0, l), &A(ui*mr, l), sizeof(T) * leni);
    }

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
                // Selected packed A panels.
                T *pakA = pakAbase + pakAsz * ui  // <packed 1st indices.
                                   + pakApn * ul; // <packed 2nd indices.
                // Starting indices.
                unsigned ist = ui*mr;
                unsigned jst = uj*nr;
                // TODO: prefetching.
                if (mker_available<T>() && leni == mr && lenj == nr)
                    ugemmn(lenk, &alpha, pakA, pakB, &beta_, &C(ist, jst), ldC);
                else if (false) // (extker_available<T>())
                    ugemmext(leni, lenj, lenk, &alpha, pakA, mr, pakB, nr, &beta_, &C(ist, jst), ldC);
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

