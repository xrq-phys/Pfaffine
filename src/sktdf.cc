/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#include "sktdf.tcc"

#ifdef PASTE_DEF
#undef PASTE_DEF
#endif
#define PASTE_DEF( funcname, typename ) \
    template signed funcname<typename> \
        (uplo_t, dim_t, typename *, inc_t, typename *, inc_t, signed *, typename *, dim_t);

PASTE_DEF( sktdf, float    )
PASTE_DEF( sktdf, double   )
PASTE_DEF( sktdf, ccscmplx )
PASTE_DEF( sktdf, ccdcmplx )
