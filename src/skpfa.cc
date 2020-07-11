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
template ccscmplx skpfa<ccscmplx>
  (char, unsigned, ccscmplx *, unsigned, unsigned,
   ccscmplx *, ccscmplx *, ccscmplx *, ccscmplx *, ccscmplx *, unsigned);
template ccdcmplx skpfa<ccdcmplx>
  (char, unsigned, ccdcmplx *, unsigned, unsigned,
   ccdcmplx *, ccdcmplx *, ccdcmplx *, ccdcmplx *, ccdcmplx *, unsigned);

