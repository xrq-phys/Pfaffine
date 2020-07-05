/**
 * \file colmaj.tcc
 * For column-major access of matrices.
 * TODO: Do margin checking;
 * TODO: Implement multi-dim.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

template <typename T>
struct colmaj
{
    T *dat;
    unsigned ld;

    colmaj() = delete;
    colmaj(T *dat_, unsigned ld_)
    : dat(dat_), ld(ld_) { }

    T &operator()(unsigned i, unsigned j)
    {
        return dat[ i + j * ld ];
    }
};
