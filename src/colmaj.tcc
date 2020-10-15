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
#pragma once
template <typename T>
struct colmaj
{
    T *dat;
    dim_t ld;

    colmaj() = delete;
    colmaj(T *dat_, dim_t ld_)
    : dat(dat_), ld(ld_) { }

    T &operator()(dim_t i, dim_t j)
    {
        return dat[ i + j * ld ];
    }
    T &operator()()
    {
        return dat[ 0 ];
    }
};
