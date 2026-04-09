.align 6
.global __prop32
.global ___prop32
.global __div65
.global ___div65
.global __mulvec
.global ___mulvec

__prop32:
___prop32:
    ldp	d1, d6, [x0, 8]
	ldr	d2, [x0]
	ldp	d4, d5, [x0, 24]
	uxtl	v2.2d, v2.2s
	ldr	d16, [x0, 120]
	uxtl	v1.2d, v1.2s
	ldp	d3, d27, [x0, 40]
	ushr	v2.2d, v2.2d, 28
	ushr	v1.2d, v1.2d, 28
	ldp	d24, d21, [x0, 56]
	uaddw	v2.2d, v2.2d, v6.2s
	uaddw	v1.2d, v1.2d, v4.2s
	mov  w7, #0xFFFF
    movk w7, #0x0FFF, LSL #16
    dup v0.4s, w7
	ushr	v6.2d, v2.2d, 28
	ushr	v4.2d, v1.2d, 28
	ldp	d20, d23, [x0, 72]
	uaddw	v6.2d, v6.2d, v5.2s
	uaddw	v4.2d, v4.2d, v3.2s
	and	v2.16b, v2.16b, v0.16b
	and	v1.16b, v1.16b, v0.16b
	ushr	v5.2d, v6.2d, 28
	ushr	v3.2d, v4.2d, 28
	and	v26.16b, v6.16b, v0.16b
	and	v25.16b, v4.16b, v0.16b
	uaddw	v5.2d, v5.2d, v27.2s
	uaddw	v3.2d, v3.2d, v24.2s
	xtn	v2.2s, v2.2d
	xtn	v1.2s, v1.2d
	ushr	v6.2d, v5.2d, 28
	ushr	v4.2d, v3.2d, 28
	ldr	q7, [x0]
	stp	d2, d1, [x0, 16]
	uaddw	v6.2d, v6.2d, v21.2s
	uaddw	v4.2d, v4.2d, v20.2s
	ldp	d22, d19, [x0, 88]
	and	v24.16b, v5.16b, v0.16b
	ushr	v2.2d, v4.2d, 28
	ushr	v5.2d, v6.2d, 28
	bic	v7.4s, #240, lsl #24
	ldp	d18, d17, [x0, 104]
	uaddw	v1.2d, v2.2d, v22.2s
	uaddw	v5.2d, v5.2d, v23.2s
	str	q7, [x0]
	and	v20.16b, v4.16b, v0.16b
	and	v7.16b, v3.16b, v0.16b
	ushr	v4.2d, v5.2d, 28
	ushr	v3.2d, v1.2d, 28
	xtn	v22.2s, v26.2d
	xtn	v2.2s, v25.2d
	uaddw	v4.2d, v4.2d, v19.2s
	uaddw	v3.2d, v3.2d, v18.2s
	and	v21.16b, v6.16b, v0.16b
	and	v6.16b, v5.16b, v0.16b
	stp	d22, d2, [x0, 32]
	and	v5.16b, v1.16b, v0.16b
	ushr	v2.2d, v4.2d, 28
	ushr	v1.2d, v3.2d, 28
	xtn	v18.2s, v24.2d
	xtn	v7.2s, v7.2d
	uaddw	v2.2d, v2.2d, v17.2s
	uaddw	v1.2d, v1.2d, v16.2s
	and	v3.16b, v3.16b, v0.16b
	and	v4.16b, v4.16b, v0.16b
	stp	d18, d7, [x0, 48]
	and	v16.16b, v2.16b, v0.16b
	and	v7.16b, v1.16b, v0.16b
	xtn	v6.2s, v6.2d
	xtn	v18.2s, v21.2d
	xtn	v17.2s, v20.2d
	xtn	v5.2s, v5.2d
	xtn	v3.2s, v3.2d
	str	d6, [x0, 80]
	xtn	v4.2s, v4.2d
	xtn	v6.2s, v7.2d
	stp	d18, d17, [x0, 64]
	xtn	v16.2s, v16.2d
	str	d5, [x0, 88]
	str	d4, [x0, 96]
	ushr	v2.2d, v2.2d, 28
	str	d3, [x0, 104]
	ushr	v1.2d, v1.2d, 28
	stp	d16, d6, [x0, 112]
	ldr	d5, [x0, 144]
	ldp	d4, d3, [x0, 128]
	ldp	d23, d22, [x0, 160]
	uaddw	v2.2d, v2.2d, v4.2s
	ldr	d4, [x0, 152]
	uaddw	v1.2d, v1.2d, v3.2s
	ldp	d21, d20, [x0, 176]
	ushr	v6.2d, v2.2d, 28
	ushr	v3.2d, v1.2d, 28
	and	v2.16b, v2.16b, v0.16b
	and	v1.16b, v1.16b, v0.16b
	uaddw	v6.2d, v6.2d, v5.2s
	uaddw	v3.2d, v3.2d, v4.2s
	xtn	v2.2s, v2.2d
	xtn	v1.2s, v1.2d
	ushr	v5.2d, v6.2d, 28
	ushr	v4.2d, v3.2d, 28
	and	v7.16b, v6.16b, v0.16b
	and	v6.16b, v3.16b, v0.16b
	stp	d2, d1, [x0, 128]
	uaddw	v5.2d, v5.2d, v23.2s
	uaddw	v4.2d, v4.2d, v22.2s
	xtn	v7.2s, v7.2d
	ldp	d19, d18, [x0, 192]
	ushr	v2.2d, v4.2d, 28
	ushr	v3.2d, v5.2d, 28
	str	d7, [x0, 144]
	and	v5.16b, v5.16b, v0.16b
	and	v4.16b, v4.16b, v0.16b
	uaddw	v1.2d, v2.2d, v20.2s
	uaddw	v3.2d, v3.2d, v21.2s
	xtn	v2.2s, v6.2d
	xtn	v5.2s, v5.2d
	ushr	v6.2d, v1.2d, 28
	ushr	v7.2d, v3.2d, 28
	and	v3.16b, v3.16b, v0.16b
	str	d2, [x0, 152]
	str	d5, [x0, 160]
	uaddw	v2.2d, v6.2d, v18.2s
	uaddw	v5.2d, v7.2d, v19.2s
	ldp	d17, d16, [x0, 208]
	and	v1.16b, v1.16b, v0.16b
	and	v6.16b, v5.16b, v0.16b
	xtn	v4.2s, v4.2d
	xtn	v3.2s, v3.2d
	and	v0.16b, v2.16b, v0.16b
	shrn	v5.2s, v5.2d, 28
	shrn	v2.2s, v2.2d, 28
	xtn	v1.2s, v1.2d
	stp	d4, d3, [x0, 168]
	xtn	v0.2s, v0.2d
	xtn	v3.2s, v6.2d
	add	v17.2s, v17.2s, v5.2s
	add	v16.2s, v16.2s, v2.2s
	str	d1, [x0, 184]
	stp	d3, d0, [x0, 192]
	stp	d17, d16, [x0, 208]
	ret

__div65:
___div65:
    ldr	q4, [x1]
	movi	v3.2s, 0x41
    mov	w2, 1009
    dup	v1.4s, w2
	ushr	v0.4s, v4.4s, 12
	movi	v2.4s, 0xf, msl 8
	mul	v1.4s, v0.4s, v1.4s
	ushr	v1.4s, v1.4s, 16
	mls	v0.4s, v1.4s, v3.s[0]
	str	q1, [x0]
	shl	v0.4s, v0.4s, 12
	and	v4.16b, v4.16b, v2.16b
	add	v0.4s, v4.4s, v0.4s
	str	q0, [x1]
	ret


___mulvec:
____mulvec:
    
    ret