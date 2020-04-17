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
template scomplex skpfa<scomplex>(char, unsigned, scomplex *, unsigned, unsigned);
template dcomplex skpfa<dcomplex>(char, unsigned, dcomplex *, unsigned, unsigned);

void sskpfa(float    *Pfa, char uplo, unsigned n, float    *A, unsigned ldA, unsigned inv)
{ *Pfa = skpfa<float>   (uplo, n, A, ldA, inv); }
void dskpfa(double   *Pfa, char uplo, unsigned n, double   *A, unsigned ldA, unsigned inv)
{ *Pfa = skpfa<double>  (uplo, n, A, ldA, inv); }
void cskpfa(scomplex *Pfa, char uplo, unsigned n, scomplex *A, unsigned ldA, unsigned inv)
{ *Pfa = skpfa<scomplex>(uplo, n, A, ldA, inv); }
void zskpfa(dcomplex *Pfa, char uplo, unsigned n, dcomplex *A, unsigned ldA, unsigned inv)
{ *Pfa = skpfa<dcomplex>(uplo, n, A, ldA, inv); }

void cal_sskpfa_(float    *Pfa, char *uplo, unsigned *n, float    *A, unsigned *ldA, unsigned *inv, 
        float    *Sp1, float    *Sp2, float    *Sp3, float    *Sp4, float    *Sp5, unsigned *npanel)
{ *Pfa = skpfa<float>   (*uplo, *n, A, *ldA, *inv, Sp1, Sp2, Sp3, Sp4, Sp5, *npanel); }
void cal_dskpfa_(double   *Pfa, char *uplo, unsigned *n, double   *A, unsigned *ldA, unsigned *inv, 
        double   *Sp1, double   *Sp2, double   *Sp3, double   *Sp4, double   *Sp5, unsigned *npanel)
{ *Pfa = skpfa<double>  (*uplo, *n, A, *ldA, *inv, Sp1, Sp2, Sp3, Sp4, Sp5, *npanel); }
void cal_cskpfa_(scomplex *Pfa, char *uplo, unsigned *n, scomplex *A, unsigned *ldA, unsigned *inv, 
        scomplex *Sp1, scomplex *Sp2, scomplex *Sp3, scomplex *Sp4, scomplex *Sp5, unsigned *npanel)
{ *Pfa = skpfa<scomplex>(*uplo, *n, A, *ldA, *inv, Sp1, Sp2, Sp3, Sp4, Sp5, *npanel); }
void cal_zskpfa_(dcomplex *Pfa, char *uplo, unsigned *n, dcomplex *A, unsigned *ldA, unsigned *inv, 
        dcomplex *Sp1, dcomplex *Sp2, dcomplex *Sp3, dcomplex *Sp4, dcomplex *Sp5, unsigned *npanel)
{ *Pfa = skpfa<dcomplex>(*uplo, *n, A, *ldA, *inv, Sp1, Sp2, Sp3, Sp4, Sp5, *npanel); }

