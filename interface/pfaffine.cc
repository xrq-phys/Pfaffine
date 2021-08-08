/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#include "pfaffine.h"
#include "pfaffine.hh"
#include "pfaffine.tcc"

#ifdef EXPAND_NAME
#undef EXPAND_NAME
#endif

// Template instantiation.
#ifdef PASTE_DEF
#undef PASTE_DEF
#endif
#define PASTE_DEF( typename ) \
    template signed skpfa< typename > \
        (char, unsigned, typename *, signed, unsigned, typename *);

PASTE_DEF( float    )
PASTE_DEF( double   )
PASTE_DEF( ccscmplx )
PASTE_DEF( ccdcmplx )
#undef EXPAND_NAME
#undef PASTE_DEF

// C simplified interface.
#define EXPAND_NAME( typechar, funcname ) typechar##funcname##_
#define PASTE_DEF( typename, typechar ) \
    void EXPAND_NAME( typechar, skpfa ) \
        (char uplo, unsigned n, typename *A, signed ldA, unsigned inv, typename *dPfa, signed *info) \
    { \
        *info = skpfa<typename>(uplo, n, A, ldA, inv, dPfa); \
    }

PASTE_DEF( float,    s )
PASTE_DEF( double,   d )
PASTE_DEF( ccscmplx, c )
PASTE_DEF( ccdcmplx, z )
#undef EXPAND_NAME
#undef PASTE_DEF

// C / Fortran full interface.
#define EXPAND_NAME( typechar, funcname ) cal_##typechar##funcname##_
#define PASTE_DEF( typename, typechar ) \
    void EXPAND_NAME( typechar, skpfa ) \
        (char *uplo,  \
         unsigned *n, \
         typename *_A, signed *ldA, \
         typename *_G, signed *ldG, \
         signed *iPiv,   \
         signed *inv,    \
         typename *dPfa, signed *info, \
         typename *_Work, unsigned *lWork) \
    { \
        *info = skpfa<typename>(check_uplo(*uplo), *n, _A, *ldA, _G, *ldG, iPiv, *inv, dPfa, _Work, *lWork); \
    }

PASTE_DEF( float,    s )
PASTE_DEF( double,   d )
PASTE_DEF( ccscmplx, c )
PASTE_DEF( ccdcmplx, z )
#undef EXPAND_NAME
#undef PASTE_DEF
