/**
 * \file floorsqrt.tcc
 * Floor of square root of an integer,
 *   credit: https://www.geeksforgeeks.org/square-root-of-an-integer.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

template <typename TI>
TI floorsqrt(TI x)
{
    // Base cases
    if (x == 0 || x == 1)
       return x;

    // Do Binary Search for floor(sqrt(x))
    int start = 1,
        end = x,
        ans;
    while (start <= end)
    {
        int mid = (start + end) / 2;

        // If x is a perfect square
        if (mid*mid == x)
            return mid;

        // Since we need floor, we update answer when mid*mid is
        // smaller than x, and move closer to sqrt(x)
        if (mid*mid < x)
        {
            start = mid + 1;
            ans = mid;
        }
        else
            // If mid*mid is greater than x
            end = mid-1;
    }
    return ans;
}

// Find i and j in flatten index of lower triangular matrix.
// NOTE: The matrix in again column-major.
template <typename TI>
void lotrg2ij(TI n, TI *i, TI *j)
{
    // n = (j - 1) * j / 2 + i;
    *j = (floorsqrt(1 + 8*n) + 1) / 2;
    *i = n - (*j) * (*j - 1) / 2;
}
