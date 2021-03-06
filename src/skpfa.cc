/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#include "skpfa.tcc"

#ifdef PASTE_DEF
#undef PASTE_DEF
#endif
#define PASTE_DEF( funcname, typename ) \
    template signed funcname<typename> \
        (uplo_t, dim_t, typename *, inc_t, typename *, inc_t, signed *, bool, typename *, typename *, dim_t);

PASTE_DEF( skpfa, float    )
PASTE_DEF( skpfa, double   )
PASTE_DEF( skpfa, ccscmplx )
PASTE_DEF( skpfa, ccdcmplx )
