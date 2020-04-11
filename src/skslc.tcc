/**
 * \file skslc.tcc
 * Get a slice from antisymmetric matrix.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#define A(i,j) A[ (i) + (j)*(ldA) ]

template<typename T>
inline void skslc(unsigned n, unsigned i, T *x, T *A, unsigned ldA)
{
    x[i] = 0.0;
    for (unsigned j = 0; j < i; ++j)
        x[j] =  A(j, i);
    for (unsigned j = i+1; j < n; ++j)
        x[j] = -A(i, j);
}

