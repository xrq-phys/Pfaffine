/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#pragma once
#include "blis.h"

inline dim_t isqr(dim_t m) { return m * m; }
inline dim_t icub(dim_t m) { return m * m * m; }

inline dim_t optpanel(dim_t n, dim_t stride)
{
    dim_t npanel = stride;

    // Perform linear search (linear is enough in most cases).
    while (2*icub(npanel) - isqr(npanel) - n/6 < 0)
        npanel += stride;

    return npanel;
}

