/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#include "skpfa.tcc"

template float skpfa<float>
  (char, unsigned, float *, unsigned, unsigned, 
   float *, float *, float *, float *, float *, unsigned);
template double skpfa<double>
  (char, unsigned, double *, unsigned, unsigned, 
   double *, double *, double *, double *, double *, unsigned);
template scomplex skpfa<scomplex>
  (char, unsigned, scomplex *, unsigned, unsigned,
   scomplex *, scomplex *, scomplex *, scomplex *, scomplex *, unsigned);
template dcomplex skpfa<dcomplex>
  (char, unsigned, dcomplex *, unsigned, unsigned,
   dcomplex *, dcomplex *, dcomplex *, dcomplex *, dcomplex *, unsigned);

