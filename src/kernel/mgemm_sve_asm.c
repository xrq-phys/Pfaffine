/*
 * \file mgemm_sve_asm.c
 * SVE kernel interface for Post-K computer (Fugaku).
 * Wrap ASM files inline.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

// Common inline ASM options.
// For security, clobber both V and Z.

#define FP_CLOBBER_SAVED \
  "v8", "v9", "v10", "v11", "v12", "v13", "v14", "v15", \
  "z8", "z9", "z10", "z11", "z12", "z13", "z14", "z15"

#define FP_CLOBBER_SP \
  "v0", "v1", "v2", "v3", "v4", "v5", "v6", "v7", \
  "v16", "v17", "v18", "v19", "v20", "v21", "v22", "v23", \
  "v24", "v25", "v26", "v27", "v28", "v29", "v30", "v31", \
  "z0", "z1", "z2", "z3", "z4", "z5", "z6", "z7", \
  "z16", "z17", "z18", "z19", "z20", "z21", "z22", "z23", \
  "z24", "z25", "z26", "z27", "z28", "z29", "z30", "z31"

#define I_CLOBBER_SAVED \
  "x16", "x17", "x18", "x19", "x20", "x21", "x22", "x23"

#define I_CLOBBER_SP \
  "x0", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", \
  "x9", "x10", "x11", "x12", "x13", "x14", "x15"

#define I_PARAMS \
  [shape]  "m" (shape),  \
  [coeffs] "m" (coeffs), \
  [A]   "m" (A),   \
  [ldA] "m" (ldA), \
  [B]   "m" (B),   \
  [ldB] "m" (ldB), \
  [C]   "m" (C),   \
  [ldC] "m" (ldC), \
  [csC] "m" (csC), \
  [nextA] "m" (nextA), \
  [nextB] "m" (nextB)

#define I_PARAMS_ASM \
  "	ldr	x0, %[shape]	\r\n" \
  "	ldr	x1, %[coeffs]	\r\n" \
  "	ldr	x2, %[A]	\r\n" \
  "	ldr	x3, %[ldA]	\r\n" \
  "	ldr	x4, %[B]	\r\n" \
  "	ldr	x5, %[ldB]	\r\n" \
  "	ldr	x6, %[C]	\r\n" \
  "	ldr	x7, %[ldC]	\r\n" \
  "	ldr	x20, %[csC]	\r\n" \
  "	ldr	x18, %[nextA]	\r\n" \
  "	ldr	x19, %[nextB]	\r\n"

// Kernels by including external ASM.

unsigned long dmgemm_2wx14(unsigned long *shape,
                           void *coeffs,
                           void *A, unsigned long ldA,
                           void *B, unsigned long ldB,
                           void *C, unsigned long ldC, long csC,
                           void *nextA, void *nextB)
{ __asm__ __volatile__ (
#ifndef _Asm_Allow_Convension
  I_PARAMS_ASM
#endif
"	.include \"kernel/dmgemm_sve2wx14.s\"\r\n"
:
:
#ifndef _Asm_Allow_Convension
  I_PARAMS
#endif
:
  FP_CLOBBER_SAVED, I_CLOBBER_SAVED
#ifndef _Asm_Allow_Convension
, FP_CLOBBER_SP, I_CLOBBER_SP
#endif
); }

unsigned long dmgemm_1wx28(unsigned long *shape,
                           void *coeffs,
                           void *A, unsigned long ldA,
                           void *B, unsigned long ldB,
                           void *C, unsigned long ldC, long csC,
                           void *nextA, void *nextB)
{ __asm__ __volatile__ (
#ifndef _Asm_Allow_Convension
  I_PARAMS_ASM
#endif
"	.include \"kernel/dmgemm_sve1wx28.s\"\r\n"
:
:
#ifndef _Asm_Allow_Convension
  I_PARAMS
#endif
:
  FP_CLOBBER_SAVED, I_CLOBBER_SAVED
#ifndef _Asm_Allow_Convension
, FP_CLOBBER_SP, I_CLOBBER_SP
#endif
); }

unsigned long zmgemm_3wx9 (unsigned long *shape,
                           double *coeffs,
                           void *A, unsigned long ldA,
                           void *B, unsigned long ldB,
                           void *C, unsigned long ldC, long csC,
                           void *nextA, void *nextB)
{ if (csC != 1) return -1; // Column striding for complex is not yet implemented.
  __asm__ __volatile__ (
#ifndef _Asm_Allow_Convension
  I_PARAMS_ASM
#endif
"	.include \"kernel/zmgemm_sve1wx28.s\"\r\n"
:
:
#ifndef _Asm_Allow_Convension
  I_PARAMS
#endif
:
  FP_CLOBBER_SAVED, I_CLOBBER_SAVED
#ifndef _Asm_Allow_Convension
, FP_CLOBBER_SP, I_CLOBBER_SP
#endif
); }

unsigned long zmgemm_2wx14(unsigned long *shape,
                           double *coeffs,
                           void *A, unsigned long ldA,
                           void *B, unsigned long ldB,
                           void *C, unsigned long ldC, long csC,
                           void *nextA, void *nextB)
{ if (csC != 1) return -1; // Column striding for complex is not yet implemented.
  __asm__ __volatile__ (
#ifndef _Asm_Allow_Convension
  I_PARAMS_ASM
#endif
"	.include \"kernel/zmgemm_sve2wx14.s\"\r\n"
:
:
#ifndef _Asm_Allow_Convension
  I_PARAMS
#endif
:
  FP_CLOBBER_SAVED, I_CLOBBER_SAVED
#ifndef _Asm_Allow_Convension
, FP_CLOBBER_SP, I_CLOBBER_SP
#endif
); }

unsigned long zmgemm_1wx28(unsigned long *shape,
                           double *coeffs,
                           void *A, unsigned long ldA,
                           void *B, unsigned long ldB,
                           void *C, unsigned long ldC, long csC,
                           void *nextA, void *nextB)
{ if (csC != 1) return -1; // Column striding for complex is not yet implemented.
  __asm__ __volatile__ (
#ifndef _Asm_Allow_Convension
  I_PARAMS_ASM
#endif
"	.include \"kernel/zmgemm_sve3wx9.s\"\r\n"
:
:
#ifndef _Asm_Allow_Convension
  I_PARAMS
#endif
:
  FP_CLOBBER_SAVED, I_CLOBBER_SAVED
#ifndef _Asm_Allow_Convension
, FP_CLOBBER_SP, I_CLOBBER_SP
#endif
); }
