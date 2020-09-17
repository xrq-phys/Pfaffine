/*
 * \file pfapack.cc
 * Interface compatible with the old Pfapack by Wimmer.
 * This file offers no header. Symbols m_?skpfa_ are injected to libpfaffine for direct calling.
 * Differences:
 *   info acts as well as input parameter, indicating whether inv is required (inv = !info).
 *   procedure names are different.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#include <iostream>
#include <cstdlib>
#include <complex>
#include "skpfa.hh"
#include "blalink.hh"
#include "pfaffine.tcc"

const static unsigned npanel = 8;

template <typename T>
void set_sp_size(unsigned n, unsigned *lWork, T *lWorkOut, signed *info)
{
    // Override also lWork.
    *lWork = npanel * n + n * n;
    // Should lie in exact precision reange.
    *lWorkOut = *lWork;
    // Exit OK.
    *info = 0;
}

signed check_sp_size(unsigned n, int nWork, int info)
{
    using namespace std;

    if (nWork < n * npanel + n * n) {
        cerr << "SKPFA: Scratchpad memory is too small." << endl;
        return err_info(Pfaffine_BAD_SCRATCHPAD, 9);
    }
    return 0;
}

template <typename T>
void la_skpfa(char uplo, char mthd, unsigned n, T *A, signed ldA, T *Pfa,
              signed *iWork, T *work, unsigned *lWork, signed *info)
{
    using namespace std;

    // Spm query.
    if (*lWork < 0)
        return set_sp_size<T>(n, lWork, work, info);
    else
        check_sp_size(n, *lWork, *info);
    switch (mthd)
    {
    case 'p':
    case 'P':
        break;

    default:
        cerr << "SKPFA: Householder transformation is not implemented. Only Parlett-Reid is supported now." << endl;
        *info = err_info(Pfaffine_NOT_IMPLEMNTED, 2);
        return;
    }

    // Set memory spaces.
    T *Sp1 = work;
    T *SpG = work + n * npanel;

    // TODO: Changed code to use only n elements.
    signed *iPov = new signed[n+1];

    // Execute SKPFA.
    *info = skpfa<T>(check_uplo(uplo), n, A, ldA, SpG, n, iPov, !*info, Pfa, Sp1, n * npanel);
}

// Instantiate.
#ifdef EXPAND_NAME
#undef EXPAND_NAME
#endif
#define EXPAND_NAME( typechar, funcname ) m_##typechar##funcname##_
#ifdef PASTE_DEF
#undef PASTE_DEF
#endif
#define PASTE_DEF( typename, typechar ) \
    extern "C" void EXPAND_NAME( typechar, skpfa ) \
        (char *uplo,  \
         char *mthd,  \
         unsigned *n, \
         typename *A, signed *ldA,     \
         typename *Pfa, signed *iWork, \
         typename *work, unsigned *lWork, signed *info) \
    { \
        la_skpfa<typename>(*uplo, *mthd, *n, A, *ldA, Pfa, iWork, work, lWork, info); \
    }

PASTE_DEF( float,    s )
PASTE_DEF( double,   d )
PASTE_DEF( ccscmplx, c )
PASTE_DEF( ccdcmplx, z )
