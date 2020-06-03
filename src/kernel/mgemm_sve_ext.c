/*
 * \file mgemm_sve_ext.c
 * SVE kernel interface for Post-K computer (Fugaku).
 * As SVE is scalable, custom-size extensions are provided.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

unsigned long dmgemm_2wx14(unsigned long *shape, 
                           double *coeffs,
                           double *A, unsigned long ldA,
                           double *B, unsigned long ldB,
                           double *C, unsigned long ldC);
unsigned long dmgemm_1wx28(unsigned long *shape, 
                           double *coeffs,
                           double *A, unsigned long ldA,
                           double *B, unsigned long ldB,
                           double *C, unsigned long ldC);
unsigned dvecln(void);

unsigned dvecln_iso(void)
{
    unsigned long vecln = 0;
    if (vecln == 0) {
        vecln = 2 * dvecln() > 14 ? 14 : 2 * dvecln();
    }
    return vecln;
}

// Using only isotropic kernels at the moment.
// TODO: Implement 14x2w kernels to utilize the full power.
void udgemm_isomax(unsigned k,
                   double alpha, double *A, double *B, 
                   double beta, double *C, unsigned ldC)
{
    // Array-parameters are static.
    const unsigned ldAB = dvecln_iso();
    static unsigned long info;
    static unsigned long shape[3];
    static double coeffs[2];
    shape[0] = dvecln_iso();
    shape[1] = dvecln_iso();
    shape[2] = k;
    coeffs[0] = alpha;
    coeffs[1] = beta;
    // Not handling error output currently.
    info = dmgemm_2wx14(shape, coeffs, A, ldAB, B, ldAB, C, ldC);
}

// As isotropic kernel is used,
// N-shape and T-shape are obviously the same.
void udgemmn(unsigned k, 
             double *Alpha_, double *A, double *B, 
             double *Beta_, double *C, unsigned ldC)
{ udgemm_isomax(k, *Alpha_, A, B, *Beta_, C, ldC); }
void udgemmt(unsigned k, 
             double *Alpha_, double *A, double *B, 
             double *Beta_, double *C, unsigned ldC)
{ udgemm_isomax(k, *Alpha_, A, B, *Beta_, C, ldC); }

// Custom-size extension.
unsigned udgemmext(unsigned m, unsigned n, unsigned k,
                   double *Alpha_, double *A, unsigned ldA, double *B, unsigned ldB,
                   double *Beta_, double *C, unsigned ldC)
{
    static unsigned long info;
    static unsigned long shape[3];
    static double coeffs[2];
    shape[0] = m;
    shape[1] = n;
    shape[2] = k;
    coeffs[0] = *Alpha_;
    coeffs[1] = *Beta_;
    if (m <= dvecln())
        info = dmgemm_1wx28(shape, coeffs, A, ldA, B, ldB, C, ldC);
    else if (m <= dvecln() * 2)
        info = dmgemm_2wx14(shape, coeffs, A, ldA, B, ldB, C, ldC);
    else
        // No such kernel.
        info = 1;
    return (unsigned)info;
}

