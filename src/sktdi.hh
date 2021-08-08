/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#pragma once
#include "blis.h"

/**
 * \brief Inverts an even-rank skew-symmetric matrix from
 *        Parlett-Reid decomposition A, G.
 * TODO: (v2.0 roadmap) G to be stored in A's lower diagonal part.
 *
 * \param uplo  Only upper-active is supported at the moment.
 * \param n     Dimension of matrix A.
 * \param _A    (IN/OUT) Base address of matrix A.
 * \param ldA   Leading dimension of A.
 * \param _G    (IN) Base address of matrix G.
 * \param ldG   Leading dimension of G.
 * \param iPiv  (IN) Pivoting reordering information, of length n+1.
 * \param _Work Scratchpad memory space, at least of length lWork.
 * \param lWork Length of scratchpad memory space.
 */
template <typename T>
signed sktdi(uplo_t uplo,
             dim_t n,
             T *_A, inc_t ldA,
             T *_G, inc_t ldG,
             signed *iPiv,
             T *_Work, dim_t lWork);
