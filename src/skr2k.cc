/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#include "skr2k.tcc"

// Block size of k.
unsigned tracblk = 8;

// Block size for GEMM.
unsigned extblk = 64;

void set_sqblk_size(unsigned n)
{ extblk = n; }

void set_tracdim_blk(unsigned n)
{ tracblk = n; }

template void skr2k<float>
  (char, char, unsigned, unsigned,
   float, float *, unsigned, float *, unsigned, float, float *, unsigned, float *);
template void skr2k<double>
  (char, char, unsigned, unsigned,
   double, double *, unsigned, double *, unsigned, double, double *, unsigned, double *);
template void skr2k<scomplex>
  (char, char, unsigned, unsigned,
   scomplex, scomplex *, unsigned, scomplex *, unsigned, scomplex, scomplex *, unsigned, scomplex *);
template void skr2k<dcomplex>
  (char, char, unsigned, unsigned,
   dcomplex, dcomplex *, unsigned, dcomplex *, unsigned, dcomplex, dcomplex *, unsigned, dcomplex *);

