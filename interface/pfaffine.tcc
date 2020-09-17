/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#include "skpfa.hh"
#include "pfaffine.hh"

// Translate uplo character into BLIS enumerates.
uplo_t check_uplo(char c)
{
    uplo_t uplo;
    switch (c)
    {
    case 'u':
    case 'U':
        return BLIS_UPPER;

    case 'l':
    case 'L':
        return BLIS_LOWER;

    default:
        return BLIS_DENSE;
    }
}

// Easy interface with automatic allocation.
template <typename T>
signed skpfa(char uplo, unsigned n, T *A, signed ldA, unsigned inv, T *dPfa)
{
    using namespace std;
    unsigned npanel = 8;

    T *Sp1 = new T[n*npanel];
    T *SpG = new T[n*n];
    signed *iPov = new signed[n+1];
    signed cInfo = skpfa<T>(check_uplo(uplo), n, A, ldA, SpG, n, iPov, inv, dPfa, Sp1, n*npanel);

    delete[] iPov;
    delete[] Sp1;
    delete[] SpG;
    return cInfo;
}
