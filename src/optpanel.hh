/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#pragma once
#include "blis.h"
#include <cmath>

inline double semisqr(dim_t m) { return pow(m, 1 + 0.2); }
inline double semicub(dim_t m) { return pow(m, 1 + 0.2 + 0.255); }

inline dim_t optpanel(dim_t n, dim_t stride)
{
    dim_t npanel = stride;

    // Perform linear search (linear is enough in most cases).
    while (semicub(npanel) - semisqr(npanel) - (double)n/6 < 0.0)
        npanel += stride;

    return npanel;
}

