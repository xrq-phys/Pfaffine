/*
 * \file findmax.tcc
 * Finds maximum absolute value of array.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#include <complex>

template<typename T>
void findmax(unsigned n, unsigned *i, T *m, T *x)
{
    using namespace std;
    // (Implicit) register vars.
    unsigned ireg;
    T mreg; 
    for (unsigned j = 0; j < n; ++j)
        if (abs(x[j]) > abs(mreg)) {
            ireg = j;
            mreg = x[j];
        }
    // Write back
    *m = mreg;
    *i = ireg;
}

