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
void findmax(unsigned n, int *i, T *m, T *x, unsigned ist=0, unsigned *iPov=nullptr)
{
    using namespace std;
    // (Implicit) register vars.
    int ireg = -1;
    T mreg = 0.0;
    for (unsigned j = ist; j < n; ++j)
        if ((!iPov || iPov[j] == j) && abs(x[j]) > abs(mreg)) {
            ireg = j;
            mreg = x[j];
        }
    // Write back
    *m = mreg;
    *i = ireg;
}

