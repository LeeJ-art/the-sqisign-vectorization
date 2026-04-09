#include <inttypes.h>
#include <encoded_sizes.h>
#include <fp2.h>
#include <arm_neon.h>

/* Arithmetic modulo X^2 + 1 */

void
fp2_set_small(fp2_t *x, const digit_t val)
{
    fp_set_small(&(x->re), val);
    fp_set_zero(&(x->im));
}

void
fp2_mul_small(fp2_t *x, const fp2_t *y, uint32_t n)
{
    fp_mul_small(&x->re, &y->re, n);
    fp_mul_small(&x->im, &y->im, n);
}

void
fp2_set_one(fp2_t *x)
{
    fp_set_one(&(x->re));
    fp_set_zero(&(x->im));
}

void
fp2_set_zero(fp2_t *x)
{
    fp_set_zero(&(x->re));
    fp_set_zero(&(x->im));
}

// Is a GF(p^2) element zero?
// Returns 0xFF...FF (true) if a=0, 0 (false) otherwise
uint32_t
fp2_is_zero(const fp2_t *a)
{
    return fp_is_zero(&(a->re)) & fp_is_zero(&(a->im));
}

// Compare two GF(p^2) elements in constant time
// Returns 0xFF...FF (true) if a=b, 0 (false) otherwise
uint32_t
fp2_is_equal(const fp2_t *a, const fp2_t *b)
{
    return fp_is_equal(&(a->re), &(b->re)) & fp_is_equal(&(a->im), &(b->im));
}

// Is a GF(p^2) element one?
// Returns 0xFF...FF (true) if a=1, 0 (false) otherwise
uint32_t
fp2_is_one(const fp2_t *a)
{
    return fp_is_equal(&(a->re), &ONE) & fp_is_zero(&(a->im));
}

void
fp2_copy(fp2_t *x, const fp2_t *y)
{
    fp_copy(&(x->re), &(y->re));
    fp_copy(&(x->im), &(y->im));
}

void
fp2_add(fp2_t *x, const fp2_t *y, const fp2_t *z)
{
    fp_add(&(x->re), &(y->re), &(z->re));
    fp_add(&(x->im), &(y->im), &(z->im));
}

void
fp2_add_one(fp2_t *x, const fp2_t *y)
{
    fp_add(&x->re, &y->re, &ONE);
    fp_copy(&x->im, &y->im);
}

void
fp2_sub(fp2_t *x, const fp2_t *y, const fp2_t *z)
{
    fp_sub(&(x->re), &(y->re), &(z->re));
    fp_sub(&(x->im), &(y->im), &(z->im));
}

void
fp2_neg(fp2_t *x, const fp2_t *y)
{
    fp_neg(&(x->re), &(y->re));
    fp_neg(&(x->im), &(y->im));
}

void
fp2_mul(fp2_t *x, const fp2_t *y, const fp2_t *z) /* F_{p^2} */
{
    fp_t t0, t1; /* gf(5* 2^248 - 1) */
    /* real, image */
    fp_add(&t0, &(y->re), &(y->im));
    fp_add(&t1, &(z->re), &(z->im));
    
    fp_mul(&t0, &t0, &t1);
    fp_mul(&t1, &(y->im), &(z->im));
    fp_mul(&(x->re), &(y->re), &(z->re));
    fp_sub(&(x->im), &t0, &t1);
    fp_sub(&(x->im), &(x->im), &(x->re));
    fp_sub(&(x->re), &(x->re), &t1);
}

void
fp2_sqr(fp2_t *x, const fp2_t *y)
{
    fp_t sum, diff;

    fp_add(&sum, &(y->re), &(y->im));
    fp_sub(&diff, &(y->re), &(y->im));
    fp_mul(&(x->im), &(y->re), &(y->im));
    fp_add(&(x->im), &(x->im), &(x->im));
    fp_mul(&(x->re), &sum, &diff);
}

void
fp2_inv(fp2_t *x)
{
    fp_t t0, t1;

    fp_sqr(&t0, &(x->re));
    fp_sqr(&t1, &(x->im));
    fp_add(&t0, &t0, &t1);
    fp_inv(&t0);
    fp_mul(&(x->re), &(x->re), &t0);
    fp_mul(&(x->im), &(x->im), &t0);
    fp_neg(&(x->im), &(x->im));
}

uint32_t
fp2_is_square(const fp2_t *x)
{
    fp_t t0, t1;

    fp_sqr(&t0, &(x->re));
    fp_sqr(&t1, &(x->im));
    fp_add(&t0, &t0, &t1);

    return fp_is_square(&t0);
}

void
fp2_sqrt(fp2_t *a)
{
    fp_t x0, x1, t0, t1;

    /* From "Optimized One-Dimensional SQIsign Verification on Intel and
     * Cortex-M4" by Aardal et al: https://eprint.iacr.org/2024/1563 */

    // x0 = \delta = sqrt(a0^2 + a1^2).
    fp_sqr(&x0, &(a->re));
    fp_sqr(&x1, &(a->im));
    fp_add(&x0, &x0, &x1);
    fp_sqrt(&x0);
    // If a1 = 0, there is a risk of \delta = -a0, which makes x0 = 0 below.
    // In that case, we restore the value \delta = a0.
    fp_select(&x0, &x0, &(a->re), fp_is_zero(&(a->im)));
    // x0 = \delta + a0, t0 = 2 * x0.
    fp_add(&x0, &x0, &(a->re));
    fp_add(&t0, &x0, &x0);

    // x1 = t0^(p-3)/4
    fp_exp3div4(&x1, &t0);

    // x0 = x0 * x1, x1 = x1 * a1, t1 = (2x0)^2.
    fp_mul(&x0, &x0, &x1);
    fp_mul(&x1, &x1, &(a->im));
    fp_add(&t1, &x0, &x0);
    fp_sqr(&t1, &t1);
    // If t1 = t0, return x0 + x1*i, otherwise x1 - x0*i.
    fp_sub(&t0, &t0, &t1);
    uint32_t f = fp_is_zero(&t0);
    fp_neg(&t1, &x0);
    fp_copy(&t0, &x1);
    fp_select(&t0, &t0, &x0, f);
    fp_select(&t1, &t1, &x1, f);

    // Check if t0 is zero
    uint32_t t0_is_zero = fp_is_zero(&t0);

    // Check whether t0, t1 are odd
    // Note: we encode to ensure canonical representation
    uint8_t tmp_bytes[FP_ENCODED_BYTES];
    fp_encode(tmp_bytes, &t0);
    uint32_t t0_is_odd = -((uint32_t)tmp_bytes[0] & 1);
    fp_encode(tmp_bytes, &t1);
    uint32_t t1_is_odd = -((uint32_t)tmp_bytes[0] & 1);

    // We negate the output if:
    // t0 is odd, or
    // t0 is zero and t1 is odd
    uint32_t negate_output = t0_is_odd | (t0_is_zero & t1_is_odd);
    fp_neg(&x0, &t0);
    fp_select(&(a->re), &t0, &x0, negate_output);
    fp_neg(&x0, &t1);
    fp_select(&(a->im), &t1, &x0, negate_output);
}

uint32_t
fp2_sqrt_verify(fp2_t *a)
{
    fp2_t t0, t1;

    fp2_copy(&t0, a);
    fp2_sqrt(a);
    fp2_sqr(&t1, a);

    return (fp2_is_equal(&t0, &t1));
}

void
fp2_half(fp2_t *x, const fp2_t *y)
{
    fp_half(&(x->re), &(y->re));
    fp_half(&(x->im), &(y->im));
}

void
fp2_batched_inv(fp2_t *x, int len)
{
    fp2_t t1[len], t2[len];
    fp2_t inverse;

    // x = x0,...,xn
    // t1 = x0, x0*x1, ... ,x0 * x1 * ... * xn
    fp2_copy(&t1[0], &x[0]);
    for (int i = 1; i < len; i++) {
        fp2_mul(&t1[i], &t1[i - 1], &x[i]);
    }

    // inverse = 1/ (x0 * x1 * ... * xn)
    fp2_copy(&inverse, &t1[len - 1]);
    fp2_inv(&inverse);

    fp2_copy(&t2[0], &inverse);
    // t2 = 1/ (x0 * x1 * ... * xn), 1/ (x0 * x1 * ... * x(n-1)) , ... , 1/xO
    for (int i = 1; i < len; i++) {
        fp2_mul(&t2[i], &t2[i - 1], &x[len - i]);
    }

    fp2_copy(&x[0], &t2[len - 1]);

    for (int i = 1; i < len; i++) {
        fp2_mul(&x[i], &t1[i - 1], &t2[len - i - 1]);
    }
}

// exponentiation using square and multiply
// Warning!! Not constant time!
void
fp2_pow_vartime(fp2_t *out, const fp2_t *x, const digit_t *exp, const int size)
{
    fp2_t acc;
    digit_t bit;

    fp2_copy(&acc, x);
    fp2_set_one(out);

    // Iterate over each word of exp
    for (int j = 0; j < size; j++) {
        // Iterate over each bit of the word
        for (int i = 0; i < RADIX; i++) {
            bit = (exp[j] >> i) & 1;
            if (bit == 1) {
                fp2_mul(out, out, &acc);
            }
            fp2_sqr(&acc, &acc);
        }
    }
}

void
fp2_print(const char *name, const fp2_t *a)
{
    printf("%s0x", name);

    uint8_t buf[FP_ENCODED_BYTES];
    fp_encode(&buf, &a->re); // Encoding ensures canonical rep
    for (int i = 0; i < FP_ENCODED_BYTES; i++) {
        printf("%02x", buf[FP_ENCODED_BYTES - i - 1]);
    }

    printf(" + i*0x");

    fp_encode(&buf, &a->im);
    for (int i = 0; i < FP_ENCODED_BYTES; i++) {
        printf("%02x", buf[FP_ENCODED_BYTES - i - 1]);
    }
    printf("\n");
}

void
fp2_encode(void *dst, const fp2_t *a)
{
    uint8_t *buf = dst;
    fp_encode(buf, &(a->re));
    fp_encode(buf + FP_ENCODED_BYTES, &(a->im));
}

uint32_t
fp2_decode(fp2_t *d, const void *src)
{
    const uint8_t *buf = src;
    uint32_t re, im;

    re = fp_decode(&(d->re), buf);
    im = fp_decode(&(d->im), buf + FP_ENCODED_BYTES);
    return re & im;
}

void
fp2_select(fp2_t *d, const fp2_t *a0, const fp2_t *a1, uint32_t ctl)
{
    fp_select(&(d->re), &(a0->re), &(a1->re), ctl);
    fp_select(&(d->im), &(a0->im), &(a1->im), ctl);
}

void
fp2_cswap(fp2_t *a, fp2_t *b, uint32_t ctl)
{
    fp_cswap(&(a->re), &(b->re), ctl);
    fp_cswap(&(a->im), &(b->im), ctl);
}

// Vectorization
uint32_t fp2_is_zero_32(const uint32x4_t* p, int x){
  uint32_t a[FP2_LIMBS];
  for (int i = 0; i < FP2_LIMBS; i++) a[i] = p[i][x];
  return fp_is_zero_32(a) & fp_is_zero_32(a+FP_LIMBS);
}

void fp2_bactched_reduction(uint32x4_t *out){
    fp_bactched_reduction(out);
    fp_bactched_reduction(out+FP_LIMBS);
}

void fp2_add_batched(uint32x4_t* out, uint32x4_t *a, uint32x4_t *b){
    fp_add_batched(out, a, b);
    fp_add_batched(out+FP_LIMBS, a+FP_LIMBS, b+FP_LIMBS);
    // reduce
    //fp2_bactched_reduction(out);
}

void fp2_sub_batched(uint32x4_t* out, uint32x4_t *a, uint32x4_t *b){
    fp_sub_batched(out, a, b);
    fp_sub_batched(out+FP_LIMBS, a+FP_LIMBS, b+FP_LIMBS);
    // reduce
    //fp2_bactched_reduction(out);
}

void fp2_mul_batched(uint32x4_t *out, uint32x4_t *a, uint32x4_t *b){
    uint32x4_t tmp[FP2_LIMBS];
    uint32x4_t ain[FP2_LIMBS], bin[FP2_LIMBS];
    memmove(ain, a, sizeof(uint32x4_t)*FP2_LIMBS);
    memmove(bin, b, sizeof(uint32x4_t)*FP2_LIMBS);

    // tmp_a = a_re + a_im
    fp_add_batched(tmp, ain, ain+FP_LIMBS);
    // tmp_b = b_re + b_im
    fp_add_batched(tmp+FP_LIMBS, bin, bin+FP_LIMBS);

    // c1 = tmp_a * tmp_b
    fp_mul_batched(((uint32x2_t*)out)+FP2_LIMBS, tmp, tmp+FP_LIMBS);

    // c2 = a_im * b_im
    fp_mul_batched(((uint32x2_t*)tmp)+FP2_LIMBS, ain+FP_LIMBS, bin+FP_LIMBS);

    // c0 = a_re * b_re
    fp_mul_batched((uint32x2_t*)out, ain, bin);

    //c1 = (a0*b1) + (a1*b0) = (a0+a1)*(b0+b1) - [(a0*b0)+(a1*b1)]
    fp_add_batched(tmp, out, tmp+FP_LIMBS);
    fp_sub_batched(out+FP_LIMBS, out+FP_LIMBS, tmp);

    //c0 = (a0*b0) - (a1*b1)
    fp_sub_batched(out, out, tmp+FP_LIMBS);

    //reduce
    fp2_bactched_reduction(out);
}

void fp2_sqr_batched(uint32x4_t* out, uint32x4_t *a){
    uint32x4_t tmp[FP2_LIMBS];
    uint32x4_t ain[FP2_LIMBS];
    memmove(ain, a, sizeof(uint32x4_t)*FP2_LIMBS);

    // tmp0 = real + img
    fp_add_batched(tmp, ain, ain+FP_LIMBS);

    // tmp1 = real - img
    fp_sub_batched(tmp+FP_LIMBS, ain, ain+FP_LIMBS);

    //c0 = (a0*a0) - (a1*a1)
    fp_mul_batched((uint32x2_t*)out, tmp, tmp+FP_LIMBS);

    // img = img + img
    fp_add_batched(tmp, ain+FP_LIMBS, ain+FP_LIMBS);

    //c1 = (a0*a1) + (a1*a0) = a0 * (2*a1)
    fp_mul_batched(((uint32x2_t*)out)+FP2_LIMBS, ain, tmp);
}

void fp2_sqrt_vec_batched_2(uint32x4_t *out, uint32x4_t *in)
{    
    uint32x4_t m1[FP_LIMBS] ={0}, m2[FP_LIMBS] = {0}, a32[FP_LIMBS] = {0}, b32[FP_LIMBS] = {0}, x32[FP_LIMBS] = {0};
    /* From "Optimized One-Dimensional SQIsign Verification on Intel and
     * Cortex-M4" by Aardal et al: https://eprint.iacr.org/2024/1563 */

    // x0 = \delta = sqrt(a0^2 + a1^2)
    // fp_sqr(&x0, &(a->re));
    // fp_sqr(&p0, &(b->re));
    // fp_sqr(&x1, &(a->im));
    // fp_sqr(&p1, &(b->im));
    for (int i=0; i<FP_LIMBS; i++){
        m1[i][0] = in[i][0];
        m1[i][1] = in[i][1];
        m1[i][2] = in[i+FP_LIMBS][0];
        m1[i][3] = in[i+FP_LIMBS][1];
    }
    fp_sqr_batched((uint32x2_t*)x32, m1);

    // fp_add(&x0, &x0, &x1);
    // fp_add(&p0, &p0, &p1);
    for (int i=0; i<FP_LIMBS; i++){
        m2[i][0] = x32[i][2];
        m2[i][1] = x32[i][3];
        a32[i] = vaddq_u32(x32[i], m2[i]);
    }
    
    // fp_sqrt(&x0);
    // fp_sqrt(&p0);
    fp_sqrt_batched(a32, a32);
    
    // If a1 = 0, there is a risk of \delta = -a0, which makes x0 = 0 below.
    // In that case, we restore the value \delta = a0.
    // fp_select(&x0, &x0, &(a->re), fp_is_zero(&(a->im)));
    // fp_select(&p0, &p0, &(b->re), fp_is_zero(&(b->im)));
    for (int i=0; i<FP_LIMBS; i++){
        m2[i][0] = in[i+FP_LIMBS][0];
        m2[i][1] = in[i+FP_LIMBS][1];
        m2[i][2] = m2[i][3] = 0;
    }
    uint32x4_t zero_flag = theta_point_is_zero(m2);
    for (int i=0; i<FP_LIMBS; i++){
        a32[i] = veorq_u32(a32[i], vandq_u32(zero_flag, veorq_u32(a32[i], m1[i])));
    }

    // x0 = \delta + a0, t0 = 2 * x0.
    // fp_add(&x0, &x0, &(a->re));
    // fp_add(&p0, &p0, &(b->re));
    for (int i=0; i<FP_LIMBS; i++){
        a32[i] = vaddq_u32(a32[i], m1[i]);
    }
    fp_bactched_reduction(a32);
    
    // fp_add(&t0, &x0, &x0);
    // fp_add(&q0, &p0, &p0);
    for (int i=0; i<FP_LIMBS; i++){
        b32[i] = vaddq_u32(a32[i], a32[i]);
        b32[i][2] = b32[i][3] = 0;
    }
    fp_bactched_reduction(b32);

    // x1 = t0^(p-3)/4
    // fp_exp3div4(&x1, &t0);
    // fp_exp3div4(&p1, &q0);
    fp_exp3div4_vec(m2, b32);
    for (int i=0; i<FP_LIMBS; i++){
        m2[i][2] = m1[i][2];
        m2[i][3] = m1[i][3];
    }

    // x0 = x0 * x1, x1 = x1 * a1, t1 = (2x0)^2.
    // fp_mul(&x0, &x0, &x1);
    // fp_mul(&p0, &p0, &p1);
    // fp_mul(&x1, &x1, &(a->im));
    // fp_mul(&p1, &p1, &(b->im));
    for (int i=0; i<FP_LIMBS; i++){
        x32[i][0] = a32[i][0];
        x32[i][1] = a32[i][1];
        x32[i][2] = m2[i][0];
        x32[i][3] = m2[i][1];
    }
    fp_mul_batched((uint32x2_t *)x32, x32, m2);

    // fp_add(&t1, &x0, &x0);
    // fp_add(&q1, &p0, &p0);
    for (int i=0; i<FP_LIMBS; i++){
        a32[i] = vaddq_u32(x32[i], x32[i]);
    }
    
    // fp_sqr(&t1, &t1);
    // fp_sqr(&q1, &q1);
    fp_sqr_batched((uint32x2_t *)a32, a32);

    // If t1 = t0, return x0 + x1*i, otherwise x1 - x0*i.
    // fp_sub(&t0, &t0, &t1);
    // fp_sub(&q0, &q0, &q1);
    // fp_neg(&t1, &x0);
    // fp_neg(&q1, &p0);
    for (int i=0; i<FP_LIMBS; i++){
        a32[i][2] = x32[i][0];
        a32[i][3] = x32[i][1];
    }
    fp_sub_batched(a32, b32, a32);
    fp_bactched_reduction(a32);

    // uint32_t f = fp_is_zero(&t0);
    // uint32_t g = fp_is_zero(&q0);
    zero_flag = theta_point_is_zero(a32);
    zero_flag[2] = zero_flag[0];
    zero_flag[3] = zero_flag[1];

    // fp_copy(&t0, &x1);
    // fp_copy(&q0, &p1);
    // fp_select(&t0, &t0, &x0, f);
    // fp_select(&q0, &q0, &p0, g);
    // fp_select(&t1, &t1, &x1, f);
    // fp_select(&q1, &q1, &p1, g);
    for (int i=0; i<FP_LIMBS; i++){
        a32[i][0] = x32[i][2];
        a32[i][1] = x32[i][3];
        a32[i] = veorq_u32(a32[i], vandq_u32(zero_flag, veorq_u32(a32[i], x32[i])));
    }

    // Check if t0 is zero
    // uint32_t t0_is_zero = fp_is_zero(&t0);
    // uint32_t q0_is_zero = fp_is_zero(&q0);
    for (int i=0; i<FP_LIMBS; i++){
        m1[i] = vdupq_n_u32(0);
    }
    zero_flag = theta_point_is_zero(a32);

    // Check whether t0, t1 are odd
    // Note: we encode to ensure canonical representation
    // int8_t tmp_bytes[FP_ENCODED_BYTES];
    // fp_encode(tmp_bytes, &t0);
    // uint32_t t0_is_odd = -((uint32_t)tmp_bytes[0] & 1);
    // fp_encode(tmp_bytes, &q0);
    // uint32_t q0_is_odd = -((uint32_t)tmp_bytes[0] & 1);
    // fp_encode(tmp_bytes, &t1);
    // uint32_t t1_is_odd = -((uint32_t)tmp_bytes[0] & 1);
    // fp_encode(tmp_bytes, &q1);
    // uint32_t q1_is_odd = -((uint32_t)tmp_bytes[0] & 1);
    m1[0] = vdupq_n_u32(1);
    fp_mul_batched((uint32x2_t *)m2, a32, m1);
    m1[0] = -(vandq_u32(m2[0], m1[0]));

    // We negate the output if:
    // t0 is odd, or
    // t0 is zero and t1 is odd
    // uint32_t negate_output = t0_is_odd | (t0_is_zero & t1_is_odd);
    // uint32_t negate_output2 = q0_is_odd | (q0_is_zero & q1_is_odd);
    zero_flag[0] = (m1[0][0] | (zero_flag[0] & m1[0][2]));
    zero_flag[1] = (m1[0][1] | (zero_flag[1] & m1[0][3]));
    zero_flag[2] = zero_flag[0];
    zero_flag[3] = zero_flag[1];

    // fp_neg(&x0, &t0);
    // fp_neg(&p0, &q0);
    // fp_neg(&x1, &t1);
    // fp_neg(&p1, &q1);
    fp_neg_vec(x32, a32);

    // fp_select(&(a->re), &t0, &x0, negate_output);
    // fp_select(&(b->re), &q0, &p0, negate_output2);
    // fp_select(&(a->im), &t1, &x1, negate_output);
    // fp_select(&(b->im), &q1, &p1, negate_output2);
    for (int i=0; i<FP_LIMBS; i++){
        a32[i] = veorq_u32(a32[i], vandq_u32(zero_flag, veorq_u32(a32[i], x32[i])));
    }

    for (int i=0; i<FP_LIMBS; i++){
        out[i][0]   = a32[i][0];
        out[i+FP_LIMBS][0] = a32[i][2];
        out[i][1]   = a32[i][1];
        out[i+FP_LIMBS][1] = a32[i][3];
    }

    return;
}

void fp2_select_vec(uint32x4_t* out, uint32x4_t* in, uint32x4_t ctl_vec){
    for (int i=0; i<FP_LIMBS; i++){
        out[i] = veorq_u32(in[i], vandq_u32(ctl_vec, veorq_u32(in[i], in[i+FP_LIMBS])));
    }
}