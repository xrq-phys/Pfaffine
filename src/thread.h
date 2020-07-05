/**
 * \file thread.h
 * For shared-memory level parallelization with OpenMP.
 * In this version OpenMP is used explicitly only in skr2k.
 * (Parallel BLAS can implicitly utilize shared-memory-parallelization.)
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#ifdef _OPENMP
#if (_OPENMP >= 200805)
#include <omp.h>
#else
#error "OpenMP version < 3.0 is not supported."
#endif
#endif

#ifndef _OPENMP
inline int omp_get_max_threads(void) { return 1; }
inline int omp_get_num_threads(void) { return 1; }
inline int omp_get_thread_num (void) { return 0; }
#endif
