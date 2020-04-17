/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#include "skpfa.hh"
#include <cstdlib>

// Easy interface with automatic allocation.
template <typename T>
T skpfa(char uplo, unsigned n, T *A, unsigned ldA, unsigned inv)
{
    using namespace std;
    unsigned npanel = 8;

    T Sp1[n*npanel];
    T Sp2[n*npanel];
    if (inv) {
        T Sp3[n*n];
        T Sp4[n*n];
        T Sp5[n*npanel];
        return skpfa<T>(uplo, n, A, ldA, 1, Sp1, Sp2, Sp3, Sp4, Sp5, npanel);
    } else
        return skpfa<T>(uplo, n, A, ldA, 0, Sp1, Sp2, 0, 0, 0, npanel);
}

