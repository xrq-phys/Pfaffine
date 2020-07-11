/**
 * \file sktdi.tcc
 * For computing inverse of tri-diagonalized antisymmetric matrix.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#include "colmaj.tcc"

template <typename T>
void sktdi(unsigned n, T *_A, unsigned ldA)
{
    colmaj<T> A(_A, ldA);
    for (unsigned i = 0; i < n; ++i)
        A(i, i) = 0.0;

    for (unsigned i = 0; i < n-1; i+=2) {
        T r = 1.0;
        for (unsigned j = i+1; j < n; j+=2) {
            r /= A(j-1, j);
            A(i, j) = -r;
            A(j, i) = r; // Inverse uses both triangles, at the moment.
            if (j != n-1) {
                r *= A(j, j+1);
                A(j+1, i) = 0.0; // Clear redundant values.
                A(i, j+1) = 0.0;
            }
        }
        if (i != n-2) {
            A(i+1, i+2) = 0.0;
            A(i+2, i+1) = 0.0;
            for (unsigned j = i+2; j < n; ++j) {
                A(j, i+1) = 0.0;
                A(i+1, j) = 0.0;
            }
        }
    }
}

