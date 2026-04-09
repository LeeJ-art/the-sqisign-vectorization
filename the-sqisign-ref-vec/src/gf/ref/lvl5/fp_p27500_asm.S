.align 6
.global __prop32
.global ___prop32
.global __div27
.global ___div27
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
	ushr	v2.2d, v2.2d, 27
	ushr	v1.2d, v1.2d, 27
	ldp	d24, d21, [x0, 56]
	uaddw	v2.2d, v2.2d, v6.2s
	uaddw	v1.2d, v1.2d, v4.2s
	mov  w7, #0xFFFF
    movk w7, #0x7FF, LSL #16
    dup v0.4s, w7
	ushr	v6.2d, v2.2d, 27
	ushr	v4.2d, v1.2d, 27
	ldp	d20, d23, [x0, 72]
	uaddw	v6.2d, v6.2d, v5.2s
	uaddw	v4.2d, v4.2d, v3.2s
	and	v2.16b, v2.16b, v0.16b
	and	v1.16b, v1.16b, v0.16b
	ushr	v5.2d, v6.2d, 27
	ushr	v3.2d, v4.2d, 27
	and	v26.16b, v6.16b, v0.16b
	and	v25.16b, v4.16b, v0.16b
	uaddw	v5.2d, v5.2d, v27.2s
	uaddw	v3.2d, v3.2d, v24.2s
	xtn	v2.2s, v2.2d
	xtn	v1.2s, v1.2d
	ushr	v6.2d, v5.2d, 27
	ushr	v4.2d, v3.2d, 27
	ldr	q7, [x0]
	stp	d2, d1, [x0, 16]
	uaddw	v6.2d, v6.2d, v21.2s
	uaddw	v4.2d, v4.2d, v20.2s
	ldp	d22, d19, [x0, 88]
	and	v24.16b, v5.16b, v0.16b
	ushr	v2.2d, v4.2d, 27
	ushr	v5.2d, v6.2d, 27
	bic	v7.4s, #248, lsl #24
	ldp	d18, d17, [x0, 104]
	uaddw	v1.2d, v2.2d, v22.2s
	uaddw	v5.2d, v5.2d, v23.2s
	str	q7, [x0]
	and	v20.16b, v4.16b, v0.16b
	and	v7.16b, v3.16b, v0.16b
	ushr	v4.2d, v5.2d, 27
	ushr	v3.2d, v1.2d, 27
	xtn	v22.2s, v26.2d
	xtn	v2.2s, v25.2d
	uaddw	v4.2d, v4.2d, v19.2s
	uaddw	v3.2d, v3.2d, v18.2s
	and	v21.16b, v6.16b, v0.16b
	and	v6.16b, v5.16b, v0.16b
	stp	d22, d2, [x0, 32]
	and	v5.16b, v1.16b, v0.16b
	ushr	v2.2d, v4.2d, 27
	ushr	v1.2d, v3.2d, 27
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
	ushr	v2.2d, v2.2d, 27
	str	d3, [x0, 104]
	ushr	v1.2d, v1.2d, 27
	stp	d16, d6, [x0, 112]
	ldp	d4, d3, [x0, 128]
	ldp	d6, d5, [x0, 144]
	uaddw	v2.2d, v2.2d, v4.2s
	uaddw	v1.2d, v1.2d, v3.2s
	ldp	d18, d17, [x0, 160]
	ushr	v4.2d, v2.2d, 27
	ushr	v3.2d, v1.2d, 27
	ldp	d28, d27, [x0, 176]
	uaddw	v4.2d, v4.2d, v6.2s
	uaddw	v3.2d, v3.2d, v5.2s
	ldp	d24, d23, [x0, 192]
	ushr	v6.2d, v4.2d, 27
	ushr	v5.2d, v3.2d, 27
	and	v2.16b, v2.16b, v0.16b
	and	v1.16b, v1.16b, v0.16b
	uaddw	v6.2d, v6.2d, v18.2s
	uaddw	v5.2d, v5.2d, v17.2s
	xtn	v2.2s, v2.2d
	xtn	v1.2s, v1.2d
	ushr	v18.2d, v6.2d, 27
	ushr	v17.2d, v5.2d, 27
	and	v25.16b, v6.16b, v0.16b
	and	v26.16b, v5.16b, v0.16b
	stp	d2, d1, [x0, 128]
	uaddw	v18.2d, v18.2d, v28.2s
	uaddw	v17.2d, v17.2d, v27.2s
	ldp	d22, d21, [x0, 208]
	ushr	v6.2d, v18.2d, 27
	ushr	v5.2d, v17.2d, 27
	and	v4.16b, v4.16b, v0.16b
	and	v3.16b, v3.16b, v0.16b
	uaddw	v6.2d, v6.2d, v24.2s
	uaddw	v5.2d, v5.2d, v23.2s
	xtn	v4.2s, v4.2d
	xtn	v3.2s, v3.2d
	ushr	v2.2d, v6.2d, 27
	ushr	v1.2d, v5.2d, 27
	ldp	d20, d19, [x0, 224]
	stp	d4, d3, [x0, 144]
	uaddw	v2.2d, v2.2d, v22.2s
	uaddw	v1.2d, v1.2d, v21.2s
	ldp	d16, d7, [x0, 240]
	ushr	v4.2d, v2.2d, 27
	ushr	v3.2d, v1.2d, 27
	and	v24.16b, v6.16b, v0.16b
	and	v23.16b, v5.16b, v0.16b
	uaddw	v4.2d, v4.2d, v20.2s
	uaddw	v3.2d, v3.2d, v19.2s
	and	v6.16b, v2.16b, v0.16b
	and	v5.16b, v1.16b, v0.16b
	ushr	v2.2d, v4.2d, 27
	ushr	v1.2d, v3.2d, 27
	and	v18.16b, v18.16b, v0.16b
	and	v17.16b, v17.16b, v0.16b
	uaddw	v2.2d, v2.2d, v16.2s
	uaddw	v1.2d, v1.2d, v7.2s
	xtn	v18.2s, v18.2d
	xtn	v17.2s, v17.2d
	and	v16.16b, v2.16b, v0.16b
	and	v7.16b, v1.16b, v0.16b
	and	v4.16b, v4.16b, v0.16b
	and	v3.16b, v3.16b, v0.16b
	stp	d18, d17, [x0, 176]
	xtn	v6.2s, v6.2d
	xtn	v5.2s, v5.2d
	xtn	v25.2s, v25.2d
	xtn	v26.2s, v26.2d
	xtn	v18.2s, v24.2d
	xtn	v17.2s, v23.2d
	stp	d6, d5, [x0, 208]
	xtn	v4.2s, v4.2d
	stp	d25, d26, [x0, 160]
	xtn	v6.2s, v7.2d
	stp	d18, d17, [x0, 192]
	xtn	v3.2s, v3.2d
	xtn	v16.2s, v16.2d
	ushr	v2.2d, v2.2d, 27
	str	d6, [x0, 248]
	ushr	v1.2d, v1.2d, 27
	stp	d4, d3, [x0, 224]
	str	d16, [x0, 240]
	ldp	d4, d3, [x0, 256]
	ldp	d16, d7, [x0, 272]
	uaddw	v2.2d, v2.2d, v4.2s
	uaddw	v1.2d, v1.2d, v3.2s
	ldp	d6, d5, [x0, 288]
	ushr	v4.2d, v2.2d, 27
	ushr	v3.2d, v1.2d, 27
	and	v2.16b, v2.16b, v0.16b
	and	v1.16b, v1.16b, v0.16b
	uaddw	v4.2d, v4.2d, v16.2s
	uaddw	v3.2d, v3.2d, v7.2s
	xtn	v2.2s, v2.2d
	xtn	v1.2s, v1.2d
	and	v7.16b, v4.16b, v0.16b
	and	v0.16b, v3.16b, v0.16b
	shrn	v4.2s, v4.2d, 27
	shrn	v3.2s, v3.2d, 27
	stp	d2, d1, [x0, 256]
	xtn	v7.2s, v7.2d
	xtn	v0.2s, v0.2d
	add	v2.2s, v6.2s, v4.2s
	add	v1.2s, v5.2s, v3.2s
	stp	d7, d0, [x0, 272]
	stp	d2, d1, [x0, 288]
	ret	

__div27:
___div27:
    ldr	q2, [x1]
	mov	w2, 1214
	movi	v3.4h, 0x1B
	dup	v1.4h, w2
	ushr	v0.4s, v2.4s, 14
	movi	v4.4s, 0x3f, msl 8
	xtn	v0.4h, v0.4s
	sqdmulh	v1.4h, v0.4h, v1.h[0]
	and	v2.16b, v2.16b, v4.16b
	mls	v0.4h, v1.4h, v3.h[0]
	uxtl	v1.4s, v1.4h
	uxtl	v0.4s, v0.4h
	shl	v0.4s, v0.4s, 14
	add	v0.4s, v0.4s, v2.4s
	str	q0, [x1]
	str	q1, [x0]
	ret

___mulvec:
____mulvec:
    
    ret