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
unsigned long zmgemm_3wx9 (unsigned long *shape,
                           double *coeffs,
                           void *A, unsigned long ldA,
                           void *B, unsigned long ldB,
                           void *C, unsigned long ldC);
unsigned long zmgemm_2wx14(unsigned long *shape,
                           double *coeffs,
                           void *A, unsigned long ldA,
                           void *B, unsigned long ldB,
                           void *C, unsigned long ldC);
unsigned long zmgemm_1wx28(unsigned long *shape,
                           double *coeffs,
                           void *A, unsigned long ldA,
                           void *B, unsigned long ldB,
                           void *C, unsigned long ldC);
unsigned dvecln(void);
unsigned zvecln(void);

unsigned vecln_iso(unsigned sve_vec_ln)
{
    unsigned long vecln = 0;
    if (vecln == 0) {
        vecln = 2 * sve_vec_ln > 14 ? 14 : 2 * sve_vec_ln;
    }
    return vecln;
}

unsigned dvecln_iso(void)
{ return vecln_iso(dvecln()); }
unsigned zvecln_iso(void)
{ return vecln_iso(zvecln()); }

unsigned zgemm_sve_mr(void)
{ return (zvecln() == 4) ? 12 : zvecln_iso(); }
unsigned zgemm_sve_nr(void)
{ return (zvecln() == 4) ? 8 : zvecln_iso(); }

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

// Complex: special 12x8 and 8x12 kernels for 512-bit vectors.
void uzgemm_12x8(unsigned k,
                 double *Alpha_, void *A, void *B,
                 double *Beta_, void *C, unsigned ldC)
{
    static unsigned long info;
    static unsigned long shape[3];
    static double coeffs[4];
    if (zvecln() != 4)
        return;
    shape[0] = 12;
    shape[1] = 8;
    shape[2] = k;
    coeffs[0] = Alpha_[0]; coeffs[1] = Alpha_[1];
    coeffs[2] = Beta_[0];  coeffs[3] = Beta_[1];
    info = zmgemm_3wx9(shape, coeffs, A, 12, B, 8, C, ldC);
}
void uzgemm_8x12(unsigned k,
                 double *Alpha_, void *A, void *B,
                 double *Beta_, void *C, unsigned ldC)
{
    static unsigned long info;
    static unsigned long shape[3];
    static double coeffs[4];
    if (zvecln() != 4)
        return;
    shape[0] = 8;
    shape[1] = 12;
    shape[2] = k;
    coeffs[0] = Alpha_[0]; coeffs[1] = Alpha_[1];
    coeffs[2] = Beta_[0];  coeffs[3] = Beta_[1];
    info = zmgemm_2wx14(shape, coeffs, A, 8, B, 12, C, ldC);
}
// Uniform fallback.
void uzgemm_isomax(unsigned k,
                   double *Alpha_, void *A, void *B,
                   double *Beta_, void *C, unsigned ldC)
{
    const unsigned ldAB = zvecln_iso();
    static unsigned long info;
    static unsigned long shape[3];
    static double coeffs[4];
    shape[0] = zvecln_iso();
    shape[1] = zvecln_iso();
    shape[2] = k;
    coeffs[0] = Alpha_[0]; coeffs[1] = Alpha_[1];
    coeffs[2] = Beta_[0];  coeffs[3] = Beta_[1];
    info = zmgemm_2wx14(shape, coeffs, A, ldAB, B, ldAB, C, ldC);
}

// As isotropic kernel is used,
// N-shape and T-shape are obviously the same.
void udgemmn(unsigned k, 
             double *Alpha_, double *A, double *B, 
             double *Beta_, double *C, unsigned ldC,
             double *nextA, double *nextB)
{ udgemm_isomax(k, *Alpha_, A, B, *Beta_, C, ldC); }
void udgemmt(unsigned k, 
             double *Alpha_, double *A, double *B, 
             double *Beta_, double *C, unsigned ldC,
             double *nextA, double *nextB)
{ udgemm_isomax(k, *Alpha_, A, B, *Beta_, C, ldC); }

// Complex situation is more complicated.
void uzgemmn(unsigned k,
             double *Alpha_, void *A, void *B,
             double *Beta_, void *C, unsigned ldC,
             void *nextA, void *nextB)
{ if (zvecln() == 4)
    uzgemm_12x8(k, Alpha_, A, B, Beta_, C, ldC);
  else uzgemm_isomax(k, Alpha_, A, B, Beta_, C, ldC); }
void uzgemmt(unsigned k,
             double *Alpha_, void *A, double *B,
             double *Beta_, void *C, unsigned ldC,
             void *nextA, void *nextB)
{ if (zvecln() == 4)
    uzgemm_8x12(k, Alpha_, A, B, Beta_, C, ldC);
  else uzgemm_isomax(k, Alpha_, A, B, Beta_, C, ldC); }

// Custom-size extension.
unsigned udgemmext(unsigned m, unsigned n, unsigned k,
                   double *Alpha_, double *A, unsigned ldA, double *B, unsigned ldB,
                   double *Beta_, double *C, unsigned ldC, double *nextA, double *nextB)
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

unsigned uzgemmext(unsigned m, unsigned n, unsigned k,
                   double *Alpha_, void *A, unsigned ldA, void *B, unsigned ldB,
                   double *Beta_, void *C, unsigned ldC, double *nextA, double *nextB)
{
    static unsigned long info;
    static unsigned long shape[3];
    static double coeffs[4];
    shape[0] = m;
    shape[1] = n;
    shape[2] = k;
    coeffs[0] = Alpha_[0]; coeffs[1] = Alpha_[1];
    coeffs[2] = Beta_[0];  coeffs[3] = Beta_[1];
    if (m <= zvecln())
        info = zmgemm_1wx28(shape, coeffs, A, ldA, B, ldB, C, ldC);
    else if (m <= zvecln() * 2)
        info = zmgemm_2wx14(shape, coeffs, A, ldA, B, ldB, C, ldC);
    else if (m <= zvecln() * 3)
        info = zmgemm_3wx9(shape, coeffs, A, ldA, B, ldB, C, ldC);
    else
        // No such kernel.
        info = 1;
    return (unsigned)info;
}

