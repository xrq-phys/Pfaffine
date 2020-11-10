/**
 * \file skpfa.tcc
 * For a Parlett-Reid method to calculate Pfaffian and inverse.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#include "skpfa.hh"
#include "sktdf.hh"
#include "sktdi.hh"
#include <iostream>
#include "blalink.hh"
#include "colmaj.tcc"

template <typename T>
signed skpfa(uplo_t uplo,
             dim_t n,
             T *_A, inc_t ldA,
             T *_G, inc_t ldG,
             signed *iPov,
             bool inv,
             T *dPfa,
             T *_Work, dim_t lWork)
{
    colmaj<T> A(_A, ldA);
    sktdf<T>(uplo, n, &A(0, 0), ldA, _G, ldG, iPov, _Work, lWork);

    signed cflp = iPov[n];
    // Pfaffian
    T PfA = 1.0;
    for (unsigned i = 0; i < n-1; i+=2)
        PfA *= A(i, i+1);
    if (inv)
        sktdi<T>(uplo, n, &A(0, 0), ldA, _G, ldG, iPov, _Work, lWork);

    // Return signed PfA.
    *dPfa = (cflp % 2) ? -PfA : PfA;

    return 0;
}
