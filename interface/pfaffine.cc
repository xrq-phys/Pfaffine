/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#include "pfaffine.h"
#include "pfaffine.hh"
#include "pfaffine.tcc"

template float    skpfa<float>   (char, unsigned, float    *, unsigned, unsigned);
template double   skpfa<double>  (char, unsigned, double   *, unsigned, unsigned);
template ccscmplx skpfa<ccscmplx>(char, unsigned, ccscmplx *, unsigned, unsigned);
template ccdcmplx skpfa<ccdcmplx>(char, unsigned, ccdcmplx *, unsigned, unsigned);

void sskpfa(float    *Pfa, char uplo, unsigned n, float    *A, unsigned ldA, unsigned inv)
{ *Pfa = skpfa<float>   (uplo, n, A, ldA, inv); }
void dskpfa(double   *Pfa, char uplo, unsigned n, double   *A, unsigned ldA, unsigned inv)
{ *Pfa = skpfa<double>  (uplo, n, A, ldA, inv); }
void cskpfa(ccscmplx *Pfa, char uplo, unsigned n, ccscmplx *A, unsigned ldA, unsigned inv)
{ *Pfa = skpfa<ccscmplx>(uplo, n, A, ldA, inv); }
void zskpfa(ccdcmplx *Pfa, char uplo, unsigned n, ccdcmplx *A, unsigned ldA, unsigned inv)
{ *Pfa = skpfa<ccdcmplx>(uplo, n, A, ldA, inv); }

void cal_sskpfa_(float    *Pfa, char *uplo, unsigned *n, float    *A, unsigned *ldA, unsigned *inv,
        float    *Sp1, float    *Sp2, float    *Sp3, float    *Sp4, float    *Sp5, unsigned *npanel)
{ *Pfa = skpfa<float>   (*uplo, *n, A, *ldA, *inv, Sp1, Sp2, Sp3, Sp4, Sp5, *npanel); }
void cal_dskpfa_(double   *Pfa, char *uplo, unsigned *n, double   *A, unsigned *ldA, unsigned *inv,
        double   *Sp1, double   *Sp2, double   *Sp3, double   *Sp4, double   *Sp5, unsigned *npanel)
{ *Pfa = skpfa<double>  (*uplo, *n, A, *ldA, *inv, Sp1, Sp2, Sp3, Sp4, Sp5, *npanel); }
void cal_cskpfa_(ccscmplx *Pfa, char *uplo, unsigned *n, ccscmplx *A, unsigned *ldA, unsigned *inv,
        ccscmplx *Sp1, ccscmplx *Sp2, ccscmplx *Sp3, ccscmplx *Sp4, ccscmplx *Sp5, unsigned *npanel)
{ *Pfa = skpfa<ccscmplx>(*uplo, *n, A, *ldA, *inv, Sp1, Sp2, Sp3, Sp4, Sp5, *npanel); }
void cal_zskpfa_(ccdcmplx *Pfa, char *uplo, unsigned *n, ccdcmplx *A, unsigned *ldA, unsigned *inv,
        ccdcmplx *Sp1, ccdcmplx *Sp2, ccdcmplx *Sp3, ccdcmplx *Sp4, ccdcmplx *Sp5, unsigned *npanel)
{ *Pfa = skpfa<ccdcmplx>(*uplo, *n, A, *ldA, *inv, Sp1, Sp2, Sp3, Sp4, Sp5, *npanel); }

