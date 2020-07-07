// For computing micro gemm for 0<M<W, N<28 and K,
//   where W is {SVE vector length in words}/2.
// Parameger registers are pointers:
// x0: =>[M,N,K]
// x1: =>[ALPHA,BETA]
// x2: =>A
// x3: LDA
// x4: =>B
// x5: LDB
// x6: =>C
// x7: LDC
// July Update 1:
// x18: =>A_NEXT
// x19: =>B_NEXT
// July Update 2:
// x20: CSC
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
// NOTE:
//  This assembly file is intended to be loaded inline.
//  All dot-directives are thus commented out.
//	.arch	armv8.2-a+sve
//	.text
//	.global	dmgemm_1wx28
//	.type	dmgemm_1wx28, %function
dmgemm_1wx28__:
//	.cfi_startproc
	ldp	x9, x10, [x0], #16	// shape M and N.
	ldr	x8, [x0]	// loop parameter K.
	mov	x11, xzr
	mov	x12, #8	// double size
	ptrue	p0.d, all	// all-true.
	whilelo	p1.d, x11, x9	// predicator.
	fmov	d0, #1.0	// exact float 1.0.
	fmov	x14, d0	// use hard float in order not to conflict with sve registers.
// Register configuration:
// Z30: A column
// Z29: B elements broadcasted
// Z28: not used
// Z[0-27]: C change buffer
	fmov	z0.d, p0/m, #0.0
	fmov	z1.d, p0/m, #0.0
	fmov	z2.d, p0/m, #0.0
	fmov	z3.d, p0/m, #0.0
	fmov	z4.d, p0/m, #0.0
	fmov	z5.d, p0/m, #0.0
	fmov	z6.d, p0/m, #0.0
	fmov	z7.d, p0/m, #0.0
	fmov	z8.d, p0/m, #0.0
	fmov	z9.d, p0/m, #0.0
	fmov	z10.d, p0/m, #0.0
	fmov	z11.d, p0/m, #0.0
	fmov	z12.d, p0/m, #0.0
	fmov	z13.d, p0/m, #0.0
	fmov	z14.d, p0/m, #0.0
	fmov	z15.d, p0/m, #0.0
	fmov	z16.d, p0/m, #0.0
	fmov	z17.d, p0/m, #0.0
	fmov	z18.d, p0/m, #0.0
	fmov	z19.d, p0/m, #0.0
	fmov	z20.d, p0/m, #0.0
	fmov	z21.d, p0/m, #0.0
	fmov	z22.d, p0/m, #0.0
	fmov	z23.d, p0/m, #0.0
	fmov	z24.d, p0/m, #0.0
	fmov	z25.d, p0/m, #0.0
	fmov	z26.d, p0/m, #0.0
	fmov	z27.d, p0/m, #0.0
K_LOOP_d1wx28:
// Load columns from A.
	ld1d	z30.d, p1/z, [x2]
	madd	x2, x3, x12, x2	// move forward
	prfm	PLDL1KEEP, [x2]	// prefetch next A column
// Apply B columns.
	madd	x16, x5, x12, x4	// calculate address in advance for prefetching
	prfm	PLDL1KEEP, [x16]	// prefetch next B column
	mov	x13, x10	// counter
	ld1rd	z29.d, p0/z, [x4, #0]	// row L column 0
	fmla	z0.d, p1/m, z30.d, z29.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW_d1wx28
	ld1rd	z29.d, p0/z, [x4, #8]	// row L column 1
	fmla	z1.d, p1/m, z30.d, z29.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW_d1wx28
	ld1rd	z29.d, p0/z, [x4, #16]	// row L column 2
	fmla	z2.d, p1/m, z30.d, z29.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW_d1wx28
	ld1rd	z29.d, p0/z, [x4, #24]	// row L column 3
	fmla	z3.d, p1/m, z30.d, z29.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW_d1wx28
	ld1rd	z29.d, p0/z, [x4, #32]	// row L column 4
	fmla	z4.d, p1/m, z30.d, z29.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW_d1wx28
	ld1rd	z29.d, p0/z, [x4, #40]	// row L column 5
	fmla	z5.d, p1/m, z30.d, z29.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW_d1wx28
	ld1rd	z29.d, p0/z, [x4, #48]	// row L column 6
	fmla	z6.d, p1/m, z30.d, z29.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW_d1wx28
	ld1rd	z29.d, p0/z, [x4, #56]	// row L column 7
	fmla	z7.d, p1/m, z30.d, z29.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW_d1wx28
	ld1rd	z29.d, p0/z, [x4, #64]	// row L column 8
	fmla	z8.d, p1/m, z30.d, z29.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW_d1wx28
	ld1rd	z29.d, p0/z, [x4, #72]	// row L column 9
	fmla	z9.d, p1/m, z30.d, z29.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW_d1wx28
	ld1rd	z29.d, p0/z, [x4, #80]	// row L column 10
	fmla	z10.d, p1/m, z30.d, z29.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW_d1wx28
	ld1rd	z29.d, p0/z, [x4, #88]	// row L column 11
	fmla	z11.d, p1/m, z30.d, z29.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW_d1wx28
	ld1rd	z29.d, p0/z, [x4, #96]	// row L column 12
	fmla	z12.d, p1/m, z30.d, z29.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW_d1wx28
	ld1rd	z29.d, p0/z, [x4, #104]	// row L column 13
	fmla	z13.d, p1/m, z30.d, z29.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW_d1wx28
	ld1rd	z29.d, p0/z, [x4, #112]	// row L column 14
	fmla	z14.d, p1/m, z30.d, z29.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW_d1wx28
	ld1rd	z29.d, p0/z, [x4, #120]	// row L column 15
	fmla	z15.d, p1/m, z30.d, z29.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW_d1wx28
	ld1rd	z29.d, p0/z, [x4, #128]	// row L column 16
	fmla	z16.d, p1/m, z30.d, z29.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW_d1wx28
	ld1rd	z29.d, p0/z, [x4, #136]	// row L column 17
	fmla	z17.d, p1/m, z30.d, z29.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW_d1wx28
	ld1rd	z29.d, p0/z, [x4, #144]	// row L column 18
	fmla	z18.d, p1/m, z30.d, z29.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW_d1wx28
	ld1rd	z29.d, p0/z, [x4, #152]	// row L column 19
	fmla	z19.d, p1/m, z30.d, z29.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW_d1wx28
	ld1rd	z29.d, p0/z, [x4, #160]	// row L column 20
	fmla	z20.d, p1/m, z30.d, z29.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW_d1wx28
	ld1rd	z29.d, p0/z, [x4, #168]	// row L column 21
	fmla	z21.d, p1/m, z30.d, z29.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW_d1wx28
	ld1rd	z29.d, p0/z, [x4, #176]	// row L column 22
	fmla	z22.d, p1/m, z30.d, z29.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW_d1wx28
	ld1rd	z29.d, p0/z, [x4, #184]	// row L column 23
	fmla	z23.d, p1/m, z30.d, z29.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW_d1wx28
	ld1rd	z29.d, p0/z, [x4, #192]	// row L column 24
	fmla	z24.d, p1/m, z30.d, z29.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW_d1wx28
	ld1rd	z29.d, p0/z, [x4, #200]	// row L column 25
	fmla	z25.d, p1/m, z30.d, z29.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW_d1wx28
	ld1rd	z29.d, p0/z, [x4, #208]	// row L column 26
	fmla	z26.d, p1/m, z30.d, z29.d
	subs	x13, x13, #1
	b.eq	NEXT_ROW_d1wx28
	ld1rd	z29.d, p0/z, [x4, #216]	// row L column 27
	fmla	z27.d, p1/m, z30.d, z29.d
//	subs	x13, x13, #1
//	b	NEXT_ROW_d1wx28
NEXT_ROW_d1wx28:
	mov	x4, x16	// move forward
	subs	x8, x8, #1
	b.ne	K_LOOP_d1wx28	// next column / row.
WRITE_MEM_d1wx28:
// Override A and B buffers:
// z[30-31]: extended alpha and beta.
// z[28-29]: C memory buffer.
	ldr	x15, [x1]	// alpha, as 64-bits
	ld1rd	z30.d, p0/z, [x1]	// alpha, to the vector.
	ld1rd	z31.d, p0/z, [x1, #8]	// beta.
// Prefetch next A and B.
	prfm	PLDL2KEEP, [x18]
	prfm	PLDL2KEEP, [x19]
// (R&)Write data back to C memory.
	cmp	x14, x15
	b.eq	UNIT_ALPHA_d1wx28
// Non-unit alpha case.
// Scale all C change buffers.
	fmul	z0.d, z0.d, z30.d
	fmul	z1.d, z1.d, z30.d
	fmul	z2.d, z2.d, z30.d
	fmul	z3.d, z3.d, z30.d
	fmul	z4.d, z4.d, z30.d
	fmul	z5.d, z5.d, z30.d
	fmul	z6.d, z6.d, z30.d
	fmul	z7.d, z7.d, z30.d
	fmul	z8.d, z8.d, z30.d
	fmul	z9.d, z9.d, z30.d
	fmul	z10.d, z10.d, z30.d
	fmul	z11.d, z11.d, z30.d
	fmul	z12.d, z12.d, z30.d
	fmul	z13.d, z13.d, z30.d
	fmul	z14.d, z14.d, z30.d
	fmul	z15.d, z15.d, z30.d
	fmul	z16.d, z16.d, z30.d
	fmul	z17.d, z17.d, z30.d
	fmul	z18.d, z18.d, z30.d
	fmul	z19.d, z19.d, z30.d
	fmul	z20.d, z20.d, z30.d
	fmul	z21.d, z21.d, z30.d
	fmul	z22.d, z22.d, z30.d
	fmul	z23.d, z23.d, z30.d
	fmul	z24.d, z24.d, z30.d
	fmul	z25.d, z25.d, z30.d
	fmul	z26.d, z26.d, z30.d
	fmul	z27.d, z27.d, z30.d
// Unit alpha case.
UNIT_ALPHA_d1wx28:
	cmp	x20, #1
	b.ne	CS_CCOL_d1wx28
// Contiguous columns.
CT_CCOL_d1wx28:
//	mov	x10, x10	// x10 itself acts as counter.
	ld1d	z28.d, p1/z, [x6]	// column vector 0
	fmad	z28.d, p1/m, z31.d, z0.d
	st1d	z28.d, p1, [x6]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6]	// column vector 1
	fmad	z28.d, p1/m, z31.d, z1.d
	st1d	z28.d, p1, [x6]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6]	// column vector 2
	fmad	z28.d, p1/m, z31.d, z2.d
	st1d	z28.d, p1, [x6]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6]	// column vector 3
	fmad	z28.d, p1/m, z31.d, z3.d
	st1d	z28.d, p1, [x6]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6]	// column vector 4
	fmad	z28.d, p1/m, z31.d, z4.d
	st1d	z28.d, p1, [x6]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6]	// column vector 5
	fmad	z28.d, p1/m, z31.d, z5.d
	st1d	z28.d, p1, [x6]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6]	// column vector 6
	fmad	z28.d, p1/m, z31.d, z6.d
	st1d	z28.d, p1, [x6]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6]	// column vector 7
	fmad	z28.d, p1/m, z31.d, z7.d
	st1d	z28.d, p1, [x6]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6]	// column vector 8
	fmad	z28.d, p1/m, z31.d, z8.d
	st1d	z28.d, p1, [x6]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6]	// column vector 9
	fmad	z28.d, p1/m, z31.d, z9.d
	st1d	z28.d, p1, [x6]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6]	// column vector 10
	fmad	z28.d, p1/m, z31.d, z10.d
	st1d	z28.d, p1, [x6]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6]	// column vector 11
	fmad	z28.d, p1/m, z31.d, z11.d
	st1d	z28.d, p1, [x6]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6]	// column vector 12
	fmad	z28.d, p1/m, z31.d, z12.d
	st1d	z28.d, p1, [x6]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6]	// column vector 13
	fmad	z28.d, p1/m, z31.d, z13.d
	st1d	z28.d, p1, [x6]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6]	// column vector 14
	fmad	z28.d, p1/m, z31.d, z14.d
	st1d	z28.d, p1, [x6]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6]	// column vector 15
	fmad	z28.d, p1/m, z31.d, z15.d
	st1d	z28.d, p1, [x6]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6]	// column vector 16
	fmad	z28.d, p1/m, z31.d, z16.d
	st1d	z28.d, p1, [x6]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6]	// column vector 17
	fmad	z28.d, p1/m, z31.d, z17.d
	st1d	z28.d, p1, [x6]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6]	// column vector 18
	fmad	z28.d, p1/m, z31.d, z18.d
	st1d	z28.d, p1, [x6]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6]	// column vector 19
	fmad	z28.d, p1/m, z31.d, z19.d
	st1d	z28.d, p1, [x6]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6]	// column vector 20
	fmad	z28.d, p1/m, z31.d, z20.d
	st1d	z28.d, p1, [x6]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6]	// column vector 21
	fmad	z28.d, p1/m, z31.d, z21.d
	st1d	z28.d, p1, [x6]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6]	// column vector 22
	fmad	z28.d, p1/m, z31.d, z22.d
	st1d	z28.d, p1, [x6]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6]	// column vector 23
	fmad	z28.d, p1/m, z31.d, z23.d
	st1d	z28.d, p1, [x6]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6]	// column vector 24
	fmad	z28.d, p1/m, z31.d, z24.d
	st1d	z28.d, p1, [x6]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6]	// column vector 25
	fmad	z28.d, p1/m, z31.d, z25.d
	st1d	z28.d, p1, [x6]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6]	// column vector 26
	fmad	z28.d, p1/m, z31.d, z26.d
	st1d	z28.d, p1, [x6]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6]	// column vector 27
	fmad	z28.d, p1/m, z31.d, z27.d
	st1d	z28.d, p1, [x6]
//	subs	x10, x10, #1
//	madd	x6, x7, x12, x6
	b	END_WRITE_MEM_d1wx28
// C has column strides.
CS_CCOL_d1wx28:
	mul	x21, x20, x12	// column stride in bytes
// Generate indices.
// Z30: index for loading C columns.
	index	z30.d, xzr, x21
	ld1d	z28.d, p1/z, [x6, z30.d]	// column vector 0
	fmad	z28.d, p1/m, z31.d, z0.d
	st1d	z28.d, p1, [x6, z30.d]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6, z30.d]	// column vector 1
	fmad	z28.d, p1/m, z31.d, z1.d
	st1d	z28.d, p1, [x6, z30.d]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6, z30.d]	// column vector 2
	fmad	z28.d, p1/m, z31.d, z2.d
	st1d	z28.d, p1, [x6, z30.d]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6, z30.d]	// column vector 3
	fmad	z28.d, p1/m, z31.d, z3.d
	st1d	z28.d, p1, [x6, z30.d]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6, z30.d]	// column vector 4
	fmad	z28.d, p1/m, z31.d, z4.d
	st1d	z28.d, p1, [x6, z30.d]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6, z30.d]	// column vector 5
	fmad	z28.d, p1/m, z31.d, z5.d
	st1d	z28.d, p1, [x6, z30.d]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6, z30.d]	// column vector 6
	fmad	z28.d, p1/m, z31.d, z6.d
	st1d	z28.d, p1, [x6, z30.d]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6, z30.d]	// column vector 7
	fmad	z28.d, p1/m, z31.d, z7.d
	st1d	z28.d, p1, [x6, z30.d]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6, z30.d]	// column vector 8
	fmad	z28.d, p1/m, z31.d, z8.d
	st1d	z28.d, p1, [x6, z30.d]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6, z30.d]	// column vector 9
	fmad	z28.d, p1/m, z31.d, z9.d
	st1d	z28.d, p1, [x6, z30.d]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6, z30.d]	// column vector 10
	fmad	z28.d, p1/m, z31.d, z10.d
	st1d	z28.d, p1, [x6, z30.d]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6, z30.d]	// column vector 11
	fmad	z28.d, p1/m, z31.d, z11.d
	st1d	z28.d, p1, [x6, z30.d]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6, z30.d]	// column vector 12
	fmad	z28.d, p1/m, z31.d, z12.d
	st1d	z28.d, p1, [x6, z30.d]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6, z30.d]	// column vector 13
	fmad	z28.d, p1/m, z31.d, z13.d
	st1d	z28.d, p1, [x6, z30.d]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6, z30.d]	// column vector 14
	fmad	z28.d, p1/m, z31.d, z14.d
	st1d	z28.d, p1, [x6, z30.d]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6, z30.d]	// column vector 15
	fmad	z28.d, p1/m, z31.d, z15.d
	st1d	z28.d, p1, [x6, z30.d]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6, z30.d]	// column vector 16
	fmad	z28.d, p1/m, z31.d, z16.d
	st1d	z28.d, p1, [x6, z30.d]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6, z30.d]	// column vector 17
	fmad	z28.d, p1/m, z31.d, z17.d
	st1d	z28.d, p1, [x6, z30.d]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6, z30.d]	// column vector 18
	fmad	z28.d, p1/m, z31.d, z18.d
	st1d	z28.d, p1, [x6, z30.d]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6, z30.d]	// column vector 19
	fmad	z28.d, p1/m, z31.d, z19.d
	st1d	z28.d, p1, [x6, z30.d]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6, z30.d]	// column vector 20
	fmad	z28.d, p1/m, z31.d, z20.d
	st1d	z28.d, p1, [x6, z30.d]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6, z30.d]	// column vector 21
	fmad	z28.d, p1/m, z31.d, z21.d
	st1d	z28.d, p1, [x6, z30.d]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6, z30.d]	// column vector 22
	fmad	z28.d, p1/m, z31.d, z22.d
	st1d	z28.d, p1, [x6, z30.d]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6, z30.d]	// column vector 23
	fmad	z28.d, p1/m, z31.d, z23.d
	st1d	z28.d, p1, [x6, z30.d]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6, z30.d]	// column vector 24
	fmad	z28.d, p1/m, z31.d, z24.d
	st1d	z28.d, p1, [x6, z30.d]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6, z30.d]	// column vector 25
	fmad	z28.d, p1/m, z31.d, z25.d
	st1d	z28.d, p1, [x6, z30.d]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6, z30.d]	// column vector 26
	fmad	z28.d, p1/m, z31.d, z26.d
	st1d	z28.d, p1, [x6, z30.d]
	subs	x10, x10, #1
	madd	x6, x7, x12, x6
	b.eq	END_WRITE_MEM_d1wx28
	ld1d	z28.d, p1/z, [x6, z30.d]	// column vector 27
	fmad	z28.d, p1/m, z31.d, z27.d
	st1d	z28.d, p1, [x6, z30.d]
//	subs	x10, x10, #1
//	madd	x6, x7, x12, x6
//	b	END_WRITE_MEM_d1wx28
// End of computation.
END_WRITE_MEM_d1wx28:
	mov	x0, #0	// return normal.
	b	END_EXEC_d1wx28
END_ERROR_d1wx28:
	mov	x0, #1	// return error.
END_EXEC_d1wx28:
//	ret
//	.cfi_endproc
