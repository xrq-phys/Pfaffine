/**
 * \file gdc.tcc
 * GDC and LCM.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

template <typename TI>
TI gcd(TI m, TI n)
{
    while (m != n)
    {
        if (m > n)
            m -= n;
        else
            n -= m;
    }
    return m;
}

template <typename TI>
TI lcm(TI m, TI n)
{ return m * n / gcd<TI>(m, n); }
