/**
 * \file opack.tcc
 * Memory packing for u?gemm procedures.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
// Macros for first-index-runs-fastest.
// FIXME: implicitly referencing macros in skr2k.tcc.
// TODO: switch to colmaj.tcc.
#define M(i,j) M[ (i) + (j)*(mr) ]

template <typename T>
void opack(unsigned n, unsigned k, unsigned stride, T *A, unsigned ldA, T *Mbase, unsigned mr)
{
    unsigned mblk = n / mr;
    unsigned mblk_= n % mr;
    if (mblk_ > 0)
        ++mblk;
    else
        mblk_ = mr;

    for (unsigned ui = 0; ui < mblk; ++ui) {
        unsigned leni = (ui+1==mblk) ? mblk_ : mr;
        T *M = Mbase + stride * ui;
        for (unsigned l = 0; l < k; ++l)
            memcpy(&M(0, l), &A(ui*mr, l), sizeof(T) * leni);
    }
}
