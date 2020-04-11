/* 
 * \file mgemm_sandy.c
 * mgemm microkernels for Sandy Bridge architecture.
 * Implemented for:
 * - double
 * Not yet implemented for:
 * - single
 * - single complex
 * - double complex
 *
 * This file is from BLIS tutorial and is subject to 3-clause BSD license:
 * 
 * Copyright (C) 2018, The University of Texas at Austin
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *  - Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  - Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *  - Neither the name(s) of the copyright holder(s) nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <mmintrin.h>
#include <xmmintrin.h>  // SSE
#include <pmmintrin.h>  // SSE2
#include <emmintrin.h>  // SSE3

/* Create macros so that the matrices are stored in column-major order */

#define A(i,j) a[ (j)*lda + (i) ]
#define B(i,j) b[ (j)*ldb + (i) ]
#define C(i,j) c[ (j)*ldc + (i) ]

typedef union
{
  __m128d v;
  double d[2];
} v2df_t;

void udgemmn( unsigned k, double *alpha_,
        double *a, double *b, double *beta_, double *c, unsigned ldc );

// Note that t is only for the size not the matrix itself.
// Hence for 4x4 kernels the N and K kernels are the same.
void udgemmt( unsigned k, double *alpha_,
        double *a, double *b, double *beta_, double *c, unsigned ldc )
{ udgemmn(k, alpha_, a, b, beta_, c, ldc); }

void udgemmn( unsigned k, double *alpha_,
        double *a, double *b, double *beta_, double *c, unsigned ldc )
{
  const unsigned lda = 4,
                 ldb = 4;
  const double one = 1.0;
  const double beta = *beta_;

  v2df_t
    alpha_vreg,
    c_00_c_10_vreg,    c_01_c_11_vreg,    c_02_c_12_vreg,    c_03_c_13_vreg,
    c_20_c_30_vreg,    c_21_c_31_vreg,    c_22_c_32_vreg,    c_23_c_33_vreg,
    a_0p_a_1p_vreg,
    a_2p_a_3p_vreg,
    b_p0_vreg, b_p1_vreg, b_p2_vreg, b_p3_vreg; 

  /* expand constant alpha */
  alpha_vreg.v = _mm_loaddup_pd( (double *) alpha_ );

  c_00_c_10_vreg.v = _mm_setzero_pd();   
  c_01_c_11_vreg.v = _mm_setzero_pd();
  c_02_c_12_vreg.v = _mm_setzero_pd(); 
  c_03_c_13_vreg.v = _mm_setzero_pd(); 
  c_20_c_30_vreg.v = _mm_setzero_pd();   
  c_21_c_31_vreg.v = _mm_setzero_pd();  
  c_22_c_32_vreg.v = _mm_setzero_pd();   
  c_23_c_33_vreg.v = _mm_setzero_pd(); 

  for ( unsigned p=0; p<k; p++ ){
    a_0p_a_1p_vreg.v = _mm_load_pd( (double *) a );
    a_2p_a_3p_vreg.v = _mm_load_pd( (double *) ( a+2 ) );
    a += lda;
    _mm_prefetch( (double *) a, 1 );

    b_p0_vreg.v = _mm_loaddup_pd( (double *) b );         /* load and duplicate */
    b_p1_vreg.v = _mm_loaddup_pd( (double *) ( b+1 ) );   /* load and duplicate */
    b_p2_vreg.v = _mm_loaddup_pd( (double *) ( b+2 ) );   /* load and duplicate */
    b_p3_vreg.v = _mm_loaddup_pd( (double *) ( b+3 ) );   /* load and duplicate */
    b += ldb;
    _mm_prefetch( (double *) b, 1 );

    /* First row and second rows */
    c_00_c_10_vreg.v += a_0p_a_1p_vreg.v * b_p0_vreg.v;
    c_01_c_11_vreg.v += a_0p_a_1p_vreg.v * b_p1_vreg.v;
    c_02_c_12_vreg.v += a_0p_a_1p_vreg.v * b_p2_vreg.v;
    c_03_c_13_vreg.v += a_0p_a_1p_vreg.v * b_p3_vreg.v;

    /* Third and fourth rows */
    c_20_c_30_vreg.v += a_2p_a_3p_vreg.v * b_p0_vreg.v;
    c_21_c_31_vreg.v += a_2p_a_3p_vreg.v * b_p1_vreg.v;
    c_22_c_32_vreg.v += a_2p_a_3p_vreg.v * b_p2_vreg.v;
    c_23_c_33_vreg.v += a_2p_a_3p_vreg.v * b_p3_vreg.v;
  }
  /* alpha factor */
  c_00_c_10_vreg.v *= alpha_vreg.v;
  c_01_c_11_vreg.v *= alpha_vreg.v;
  c_02_c_12_vreg.v *= alpha_vreg.v;
  c_03_c_13_vreg.v *= alpha_vreg.v;

  c_20_c_30_vreg.v *= alpha_vreg.v;
  c_21_c_31_vreg.v *= alpha_vreg.v;
  c_22_c_32_vreg.v *= alpha_vreg.v;
  c_23_c_33_vreg.v *= alpha_vreg.v;

  /* Fetch the next matrix */
  // _mm_prefetch( (double *) anext, 1 );
  // _mm_prefetch( (double *) bnext, 1 );
  // _mm_prefetch( (double *) cnext, 2 );

  /* write back to memory */
  if (beta == one) {
    C( 0, 0 ) += c_00_c_10_vreg.d[0];  C( 0, 1 ) += c_01_c_11_vreg.d[0];
    C( 0, 2 ) += c_02_c_12_vreg.d[0];  C( 0, 3 ) += c_03_c_13_vreg.d[0];
    C( 1, 0 ) += c_00_c_10_vreg.d[1];  C( 1, 1 ) += c_01_c_11_vreg.d[1];
    C( 1, 2 ) += c_02_c_12_vreg.d[1];  C( 1, 3 ) += c_03_c_13_vreg.d[1];
    C( 2, 0 ) += c_20_c_30_vreg.d[0];  C( 2, 1 ) += c_21_c_31_vreg.d[0];
    C( 2, 2 ) += c_22_c_32_vreg.d[0];  C( 2, 3 ) += c_23_c_33_vreg.d[0];
    C( 3, 0 ) += c_20_c_30_vreg.d[1];  C( 3, 1 ) += c_21_c_31_vreg.d[1];
    C( 3, 2 ) += c_22_c_32_vreg.d[1];  C( 3, 3 ) += c_23_c_33_vreg.d[1];
  } else {
    C( 0, 0 ) = C( 0, 0 )*beta + c_00_c_10_vreg.d[0];  C( 0, 1 ) = C( 0, 1 )*beta + c_01_c_11_vreg.d[0];
    C( 0, 2 ) = C( 0, 2 )*beta + c_02_c_12_vreg.d[0];  C( 0, 3 ) = C( 0, 3 )*beta + c_03_c_13_vreg.d[0];
    C( 1, 0 ) = C( 1, 0 )*beta + c_00_c_10_vreg.d[1];  C( 1, 1 ) = C( 1, 1 )*beta + c_01_c_11_vreg.d[1];
    C( 1, 2 ) = C( 1, 2 )*beta + c_02_c_12_vreg.d[1];  C( 1, 3 ) = C( 1, 3 )*beta + c_03_c_13_vreg.d[1];
    C( 2, 0 ) = C( 2, 0 )*beta + c_20_c_30_vreg.d[0];  C( 2, 1 ) = C( 2, 1 )*beta + c_21_c_31_vreg.d[0];
    C( 2, 2 ) = C( 2, 2 )*beta + c_22_c_32_vreg.d[0];  C( 2, 3 ) = C( 2, 3 )*beta + c_23_c_33_vreg.d[0];
    C( 3, 0 ) = C( 3, 0 )*beta + c_20_c_30_vreg.d[1];  C( 3, 1 ) = C( 3, 1 )*beta + c_21_c_31_vreg.d[1];
    C( 3, 2 ) = C( 3, 2 )*beta + c_22_c_32_vreg.d[1];  C( 3, 3 ) = C( 3, 3 )*beta + c_23_c_33_vreg.d[1];
  }
}

