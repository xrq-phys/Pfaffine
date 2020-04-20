// Gets vector size for SVE.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
	.arch	armv8.2-a+sve
	.text
	.global	svecln_get 
	.type	svecln_get, %function
svecln_get:
	.cfi_startproc
	mov	x0, #0
	incw	x0
	ret
	.cfi_endproc
