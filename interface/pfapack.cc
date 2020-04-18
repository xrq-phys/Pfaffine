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

const static unsigned npanel = 8;

void set_sp_size(unsigned n, int *lWork, int *info)
{
    if (*info) {
        // Pfaffian only.
        *lWork = npanel * n * 2;
    } else {
        // With inverse.
        *lWork = npanel * n * 3 + n * n * 2;
    }
    *info = 0;
}

void check_sp_size(unsigned n, int nWork, int info)
{
    using namespace std;

    if (( info && nWork < 2 * n * npanel) || 
        (!info && nWork < 3 * n * npanel + 2 * n * n)) {
        cerr << "Pfaffine error: scratchpad memory is too small." << endl;
        _Exit(EXIT_FAILURE);
    }
}

template <typename T>
void la_skpfa(char *uplo, char *mthd, unsigned *n, T *A, unsigned *ldA, T *Pfa, 
              int *iWork, T *work, int *lWork, int *info)
{
    using namespace std;
    // Spm query.
    if (*lWork < 0)
        return set_sp_size(*n, lWork, info);
    else
        check_sp_size(*n, *lWork, *info);

    // Set memory spaces.
    T *Sp1 = work;
    T *Sp2 = work + 1 * (*n) * npanel;
    T *Sp3 = work + 2 * (*n) * npanel;
    T *Sp4 = Sp3 + 1 * (*n) * (*n);
    T *Sp5 = Sp3 + 2 * (*n) * (*n);

    // Execute
    *Pfa = skpfa<T>(*uplo, *n, A, *ldA, !(*info), Sp1, Sp2, Sp3, Sp4, Sp5, npanel);
    // TODO: Support error as code.
    *info = 0;
}

// Instantiate.
extern "C" void m_sskpfa_(char *uplo, char *mthd, unsigned *n, float *A, unsigned *ldA, float *Pfa, 
                          int *iWork, float *work, int *lWork, int *info)
{ la_skpfa<float>(uplo, mthd, n, A, ldA, Pfa, iWork, work, lWork, info); }
extern "C" void m_dskpfa_(char *uplo, char *mthd, unsigned *n, double *A, unsigned *ldA, double *Pfa, 
                          int *iWork, double *work, int *lWork, int *info)
{ la_skpfa<double>(uplo, mthd, n, A, ldA, Pfa, iWork, work, lWork, info); }

extern "C" void m_cskpfa_(char *uplo, char *mthd, unsigned *n, void *A, unsigned *ldA, void *Pfa, 
                          int *iWork, void *work, int *lWork, int *info)
{ la_skpfa<std::complex<float> >(uplo, mthd, n, (std::complex<float> *)A, ldA, (std::complex<float> *)Pfa, 
                                 iWork, (std::complex<float> *)work, lWork, info); 
}
extern "C" void m_zskpfa_(char *uplo, char *mthd, unsigned *n, void *A, unsigned *ldA, void *Pfa, 
                          int *iWork, void *work, int *lWork, int *info)
{ la_skpfa<std::complex<double> >(uplo, mthd, n, (std::complex<double> *)A, ldA, (std::complex<double> *)Pfa, 
                                  iWork, (std::complex<double> *)work, lWork, info); 
}

