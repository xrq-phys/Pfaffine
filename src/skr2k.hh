/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#pragma once

template<typename T>
void skr2k(char uplo, char trans, unsigned n, unsigned k,
           T alpha, T *A, unsigned ldA, T *B, unsigned ldB, T beta, T *C, unsigned ldC);

