/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#pragma once
#include "blis.h"

/**
 * \brief Tri-diagonal factorization of even-ranked skew-symmetric matrix A.
 * TODO: (v2.0 roadmap) It would be possible to store G's information into A's
 *       lower triangular part so that parameter G can be omitted.
 *
 * \param uplo  Only upper-active is supported at the moment.
 * \param n     Dimension of matrix A.
 * \param _A    (IN/OUT) Base address of matrix A.
 * \param ldA   Leading dimension of A.
 * \param _G    (OUT) Base address of matrix G, used to store Gaussian transformations.
 * \param ldG   Leading dimension of G.
 * \param iPov  (OUT) Pivoting reordering information.
 * \param _Work Scratchpad memory space, at least of length lWork.
 * \param lWork Length of scratchpad memory space.
 */
template <typename T>
signed sktdf(uplo_t uplo,
             dim_t n,
             T *_A, inc_t ldA,
             T *_G, inc_t ldG,
             signed *iPov,
             T *_Work, dim_t lWork);
