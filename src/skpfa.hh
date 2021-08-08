/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#pragma once
#include "blis.h"

/**
 * \brief Calculate Pfaffian (&optionally inverse) of antisymmetric matrix A.
 *        The tri-diagonal decomposition is returned after computation is complete.
 *
 * \param uplo  Only 'U' allowed.
 * \param n     Dimension of A.
 * \param _A    (IN/OUT) Array (base address of) A. Contains tri-diag part of decomposition on exit.
 * \param ldA   Leading dimension size of A (row-skip).
 * \param _G    (OUT) Array of the same size as A. Contains gaussian elimination part on exit.
 * \param ldG   Leading dimension size of G (row-skip).
 * \param iPiv  (OUT) Povoting information of size n+1.
 * \param inv   Whether inverse would be calculated and stored in A.
 * \param dPfa  (OUT) Single-FP buffer for returning final Pfaffian computed.
 * \param _Work Scratch space. Use 8*n if not sure.
 * \param lWork Scratch space size.
 */
template <typename T>
signed skpfa(uplo_t uplo,
             dim_t n,
             T *_A, inc_t ldA,
             T *_G, inc_t ldG,
             signed *iPiv,
             bool inv,
             T *dPfa,
             T *_Work, dim_t lWork);
