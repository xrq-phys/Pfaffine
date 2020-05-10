/*
 * \file vecln_sve.c
 * SVE determine vector length.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

unsigned long svecln_get(void);
unsigned long dvecln_get(void);

unsigned svecln(void)
{
    static unsigned long vecln = 0;
    if (vecln == 0)
        vecln = svecln_get();
    return vecln;
}

unsigned dvecln(void)
{
    static unsigned long vecln = 0;
    if (vecln == 0)
        vecln = dvecln_get();
    return vecln;
}

unsigned cvecln(void)
{ return svecln() / 2; }

unsigned zvecln(void)
{ return dvecln() / 2; }

