// clang-format off
// Command line : python monty.py 64
// 0x40ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

#include <stdint.h>
#include <stdio.h>

#define sspint int64_t
#define spint uint64_t
#define udpint __uint128_t
#define dpint __uint128_t

#define Wordlength 64
#define Nlimbs 7
#define Radix 55
#define Nbits 383
#define Nbytes 48

#define MONTGOMERY
// propagate carries
inline static spint prop(spint *n) {
  int i;
  spint mask = ((spint)1 << 55u) - (spint)1;
  sspint carry = (sspint)n[0];
  carry >>= 55u;
  n[0] &= mask;
  for (i = 1; i < 6; i++) {
    carry += (sspint)n[i];
    n[i] = (spint)carry & mask;
    carry >>= 55u;
  }
  n[6] += (spint)carry;
  return -((n[6] >> 1) >> 62u);
}

// propagate carries and add p if negative, propagate carries again
inline static int flatten(spint *n) {
  spint carry = prop(n);
  n[0] -= (spint)1u & carry;
  n[6] += ((spint)0x10400000000000u) & carry;
  (void)prop(n);
  return (int)(carry & 1);
}

// Montgomery final subtract
inline static int modfsb(spint *n) {
  n[0] += (spint)1u;
  n[6] -= (spint)0x10400000000000u;
  return flatten(n);
}

// Modular addition - reduce less than 2p
inline static void modadd(const spint *a, const spint *b, spint *n) {
  spint carry;
  n[0] = a[0] + b[0];
  n[1] = a[1] + b[1];
  n[2] = a[2] + b[2];
  n[3] = a[3] + b[3];
  n[4] = a[4] + b[4];
  n[5] = a[5] + b[5];
  n[6] = a[6] + b[6];
  n[0] += (spint)2u;
  n[6] -= (spint)0x20800000000000u;
  carry = prop(n);
  n[0] -= (spint)2u & carry;
  n[6] += ((spint)0x20800000000000u) & carry;
  (void)prop(n);
}

// Modular subtraction - reduce less than 2p
inline static void modsub(const spint *a, const spint *b, spint *n) {
  spint carry;
  n[0] = a[0] - b[0];
  n[1] = a[1] - b[1];
  n[2] = a[2] - b[2];
  n[3] = a[3] - b[3];
  n[4] = a[4] - b[4];
  n[5] = a[5] - b[5];
  n[6] = a[6] - b[6];
  carry = prop(n);
  n[0] -= (spint)2u & carry;
  n[6] += ((spint)0x20800000000000u) & carry;
  (void)prop(n);
}

// Modular negation
inline static void modneg(const spint *b, spint *n) {
  spint carry;
  n[0] = (spint)0 - b[0];
  n[1] = (spint)0 - b[1];
  n[2] = (spint)0 - b[2];
  n[3] = (spint)0 - b[3];
  n[4] = (spint)0 - b[4];
  n[5] = (spint)0 - b[5];
  n[6] = (spint)0 - b[6];
  carry = prop(n);
  n[0] -= (spint)2u & carry;
  n[6] += ((spint)0x20800000000000u) & carry;
  (void)prop(n);
}

// Overflow limit   = 340282366920938463463374607431768211456
// maximum possible = 9251314080475062396111552646217735
// Modular multiplication, c=a*b mod 2p
inline static void modmul(const spint *a, const spint *b, spint *c) {
  dpint t = 0;
  spint p6 = 0x10400000000000u;
  spint q = ((spint)1 << 55u); // q is unsaturated radix
  spint mask = (spint)(q - (spint)1);
  t += (dpint)a[0] * b[0];
  spint v0 = ((spint)t & mask);
  t >>= 55;
  t += (dpint)a[0] * b[1];
  t += (dpint)a[1] * b[0];
  spint v1 = ((spint)t & mask);
  t >>= 55;
  t += (dpint)a[0] * b[2];
  t += (dpint)a[1] * b[1];
  t += (dpint)a[2] * b[0];
  spint v2 = ((spint)t & mask);
  t >>= 55;
  t += (dpint)a[0] * b[3];
  t += (dpint)a[1] * b[2];
  t += (dpint)a[2] * b[1];
  t += (dpint)a[3] * b[0];
  spint v3 = ((spint)t & mask);
  t >>= 55;
  t += (dpint)a[0] * b[4];
  t += (dpint)a[1] * b[3];
  t += (dpint)a[2] * b[2];
  t += (dpint)a[3] * b[1];
  t += (dpint)a[4] * b[0];
  spint v4 = ((spint)t & mask);
  t >>= 55;
  t += (dpint)a[0] * b[5];
  t += (dpint)a[1] * b[4];
  t += (dpint)a[2] * b[3];
  t += (dpint)a[3] * b[2];
  t += (dpint)a[4] * b[1];
  t += (dpint)a[5] * b[0];
  spint v5 = ((spint)t & mask);
  t >>= 55;
  t += (dpint)a[0] * b[6];
  t += (dpint)a[1] * b[5];
  t += (dpint)a[2] * b[4];
  t += (dpint)a[3] * b[3];
  t += (dpint)a[4] * b[2];
  t += (dpint)a[5] * b[1];
  t += (dpint)a[6] * b[0];
  t += (dpint)v0 * (dpint)p6;
  spint v6 = ((spint)t & mask);
  t >>= 55;
  t += (dpint)a[1] * b[6];
  t += (dpint)a[2] * b[5];
  t += (dpint)a[3] * b[4];
  t += (dpint)a[4] * b[3];
  t += (dpint)a[5] * b[2];
  t += (dpint)a[6] * b[1];
  t += (dpint)v1 * (dpint)p6;
  c[0] = ((spint)t & mask);
  t >>= 55;
  t += (dpint)a[2] * b[6];
  t += (dpint)a[3] * b[5];
  t += (dpint)a[4] * b[4];
  t += (dpint)a[5] * b[3];
  t += (dpint)a[6] * b[2];
  t += (dpint)v2 * (dpint)p6;
  c[1] = ((spint)t & mask);
  t >>= 55;
  t += (dpint)a[3] * b[6];
  t += (dpint)a[4] * b[5];
  t += (dpint)a[5] * b[4];
  t += (dpint)a[6] * b[3];
  t += (dpint)v3 * (dpint)p6;
  c[2] = ((spint)t & mask);
  t >>= 55;
  t += (dpint)a[4] * b[6];
  t += (dpint)a[5] * b[5];
  t += (dpint)a[6] * b[4];
  t += (dpint)v4 * (dpint)p6;
  c[3] = ((spint)t & mask);
  t >>= 55;
  t += (dpint)a[5] * b[6];
  t += (dpint)a[6] * b[5];
  t += (dpint)v5 * (dpint)p6;
  c[4] = ((spint)t & mask);
  t >>= 55;
  t += (dpint)a[6] * b[6];
  t += (dpint)v6 * (dpint)p6;
  c[5] = ((spint)t & mask);
  t >>= 55;
  c[6] = (spint)t;
}

// Modular squaring, c=a*a  mod 2p
inline static void modsqr(const spint *a, spint *c) {
  udpint tot;
  udpint t = 0;
  spint p6 = 0x10400000000000u;
  spint q = ((spint)1 << 55u); // q is unsaturated radix
  spint mask = (spint)(q - (spint)1);
  tot = (udpint)a[0] * a[0];
  t = tot;
  spint v0 = ((spint)t & mask);
  t >>= 55;
  tot = (udpint)a[0] * a[1];
  tot *= 2;
  t += tot;
  spint v1 = ((spint)t & mask);
  t >>= 55;
  tot = (udpint)a[0] * a[2];
  tot *= 2;
  tot += (udpint)a[1] * a[1];
  t += tot;
  spint v2 = ((spint)t & mask);
  t >>= 55;
  tot = (udpint)a[0] * a[3];
  tot += (udpint)a[1] * a[2];
  tot *= 2;
  t += tot;
  spint v3 = ((spint)t & mask);
  t >>= 55;
  tot = (udpint)a[0] * a[4];
  tot += (udpint)a[1] * a[3];
  tot *= 2;
  tot += (udpint)a[2] * a[2];
  t += tot;
  spint v4 = ((spint)t & mask);
  t >>= 55;
  tot = (udpint)a[0] * a[5];
  tot += (udpint)a[1] * a[4];
  tot += (udpint)a[2] * a[3];
  tot *= 2;
  t += tot;
  spint v5 = ((spint)t & mask);
  t >>= 55;
  tot = (udpint)a[0] * a[6];
  tot += (udpint)a[1] * a[5];
  tot += (udpint)a[2] * a[4];
  tot *= 2;
  tot += (udpint)a[3] * a[3];
  t += tot;
  t += (udpint)v0 * p6;
  spint v6 = ((spint)t & mask);
  t >>= 55;
  tot = (udpint)a[1] * a[6];
  tot += (udpint)a[2] * a[5];
  tot += (udpint)a[3] * a[4];
  tot *= 2;
  t += tot;
  t += (udpint)v1 * p6;
  c[0] = ((spint)t & mask);
  t >>= 55;
  tot = (udpint)a[2] * a[6];
  tot += (udpint)a[3] * a[5];
  tot *= 2;
  tot += (udpint)a[4] * a[4];
  t += tot;
  t += (udpint)v2 * p6;
  c[1] = ((spint)t & mask);
  t >>= 55;
  tot = (udpint)a[3] * a[6];
  tot += (udpint)a[4] * a[5];
  tot *= 2;
  t += tot;
  t += (udpint)v3 * p6;
  c[2] = ((spint)t & mask);
  t >>= 55;
  tot = (udpint)a[4] * a[6];
  tot *= 2;
  tot += (udpint)a[5] * a[5];
  t += tot;
  t += (udpint)v4 * p6;
  c[3] = ((spint)t & mask);
  t >>= 55;
  tot = (udpint)a[5] * a[6];
  tot *= 2;
  t += tot;
  t += (udpint)v5 * p6;
  c[4] = ((spint)t & mask);
  t >>= 55;
  tot = (udpint)a[6] * a[6];
  t += tot;
  t += (udpint)v6 * p6;
  c[5] = ((spint)t & mask);
  t >>= 55;
  c[6] = (spint)t;
}

// copy
inline static void modcpy(const spint *a, spint *c) {
  int i;
  for (i = 0; i < 7; i++) {
    c[i] = a[i];
  }
}

// square n times
static void modnsqr(spint *a, int n) {
  int i;
  for (i = 0; i < n; i++) {
    modsqr(a, a);
  }
}

// Calculate progenitor
static void modpro(const spint *w, spint *z) {
  spint x[7];
  spint t0[7];
  spint t1[7];
  spint t2[7];
  spint t3[7];
  spint t4[7];
  spint t5[7];
  modcpy(w, x);
  modsqr(x, z);
  modsqr(z, t0);
  modmul(x, t0, t1);
  modmul(z, t1, z);
  modsqr(z, t0);
  modsqr(t0, t3);
  modsqr(t3, t4);
  modsqr(t4, t2);
  modcpy(t2, t5);
  modnsqr(t5, 3);
  modmul(t2, t5, t2);
  modcpy(t2, t5);
  modnsqr(t5, 6);
  modmul(t2, t5, t2);
  modcpy(t2, t5);
  modnsqr(t5, 2);
  modmul(t4, t5, t5);
  modnsqr(t5, 13);
  modmul(t2, t5, t2);
  modcpy(t2, t5);
  modnsqr(t5, 2);
  modmul(t4, t5, t4);
  modnsqr(t4, 28);
  modmul(t2, t4, t2);
  modsqr(t2, t4);
  modmul(t3, t4, t3);
  modnsqr(t3, 59);
  modmul(t2, t3, t2);
  modmul(t1, t2, t1);
  modmul(z, t1, z);
  modmul(t0, z, t0);
  modmul(t1, t0, t1);
  modsqr(t1, t2);
  modmul(t1, t2, t2);
  modsqr(t2, t2);
  modmul(t1, t2, t2);
  modmul(t0, t2, t0);
  modmul(z, t0, z);
  modsqr(z, t2);
  modmul(z, t2, t2);
  modmul(t0, t2, t0);
  modmul(t1, t0, t1);
  modcpy(t1, t2);
  modnsqr(t2, 128);
  modmul(t1, t2, t1);
  modmul(t0, t1, t0);
  modnsqr(t0, 125);
  modmul(z, t0, z);
}

// calculate inverse, provide progenitor h if available
static void modinv(const spint *x, const spint *h, spint *z) {
  spint s[7];
  spint t[7];
  if (h == NULL) {
    modpro(x, t);
  } else {
    modcpy(h, t);
  }
  modcpy(x, s);
  modnsqr(t, 2);
  modmul(s, t, z);
}

// Convert m to n-residue form, n=nres(m)
static void nres(const spint *m, spint *n) {
  const spint c[7] = {0xfc0fc0fc0fc4du,  0x781f81f81f81f8u, 0x3f03f03f03f03u,
                      0x7e07e07e07e07eu, 0x40fc0fc0fc0fc0u, 0x1f81f81f81f81fu,
                      0xcff03f03f03f0u};
  modmul(m, c, n);
}

// Convert n back to normal form, m=redc(n)
static void redc(const spint *n, spint *m) {
  int i;
  spint c[7];
  c[0] = 1;
  for (i = 1; i < 7; i++) {
    c[i] = 0;
  }
  modmul(n, c, m);
  (void)modfsb(m);
}

// is unity?
static int modis1(const spint *a) {
  int i;
  spint c[7];
  spint c0;
  spint d = 0;
  redc(a, c);
  for (i = 1; i < 7; i++) {
    d |= c[i];
  }
  c0 = (spint)c[0];
  return ((spint)1 & ((d - (spint)1) >> 55u) &
          (((c0 ^ (spint)1) - (spint)1) >> 55u));
}

// is zero?
static int modis0(const spint *a) {
  int i;
  spint c[7];
  spint d = 0;
  redc(a, c);
  for (i = 0; i < 7; i++) {
    d |= c[i];
  }
  return ((spint)1 & ((d - (spint)1) >> 55u));
}

// set to zero
static void modzer(spint *a) {
  int i;
  for (i = 0; i < 7; i++) {
    a[i] = 0;
  }
}

// set to one
static void modone(spint *a) {
  int i;
  a[0] = 1;
  for (i = 1; i < 7; i++) {
    a[i] = 0;
  }
  nres(a, a);
}

// set to integer
static void modint(int x, spint *a) {
  int i;
  a[0] = (spint)x;
  for (i = 1; i < 7; i++) {
    a[i] = 0;
  }
  nres(a, a);
}

// Modular multiplication by an integer, c=a*b mod 2p
inline static void modmli(const spint *a, int b, spint *c) {
  spint t[7];
  modint(b, t);
  modmul(a, t, c);
}

// Test for quadratic residue
static int modqr(const spint *h, const spint *x) {
  spint r[7];
  if (h == NULL) {
    modpro(x, r);
    modsqr(r, r);
  } else {
    modsqr(h, r);
  }
  modmul(r, x, r);
  return modis1(r) | modis0(x);
}

// conditional move g to f if d=1
// strongly recommend inlining be disabled using compiler specific syntax
static void modcmv(int b, const spint *g, volatile spint *f) {
  int i;
  spint c0, c1, s, t;
  spint r = 0x3cc3c33c5aa5a55au;
  c0 = (1 - b) + r;
  c1 = b + r;
  for (i = 0; i < 7; i++) {
    s = g[i];
    t = f[i];
    f[i] = c0 * t + c1 * s;
    f[i] -= r * (t + s);
  }
}

// conditional swap g and f if d=1
// strongly recommend inlining be disabled using compiler specific syntax
static void modcsw(int b, volatile spint *g, volatile spint *f) {
  int i;
  spint c0, c1, s, t, w;
  spint r = 0x3cc3c33c5aa5a55au;
  c0 = (1 - b) + r;
  c1 = b + r;
  for (i = 0; i < 7; i++) {
    s = g[i];
    t = f[i];
    w = r * (t + s);
    f[i] = c0 * t + c1 * s;
    f[i] -= w;
    g[i] = c0 * s + c1 * t;
    g[i] -= w;
  }
}

// Modular square root, provide progenitor h if available, NULL if not
static void modsqrt(const spint *x, const spint *h, spint *r) {
  spint s[7];
  spint y[7];
  if (h == NULL) {
    modpro(x, y);
  } else {
    modcpy(h, y);
  }
  modmul(y, x, s);
  modcpy(s, r);
}

// shift left by less than a word
static void modshl(unsigned int n, spint *a) {
  int i;
  a[6] = ((a[6] << n)) | (a[5] >> (55u - n));
  for (i = 5; i > 0; i--) {
    a[i] = ((a[i] << n) & (spint)0x7fffffffffffff) | (a[i - 1] >> (55u - n));
  }
  a[0] = (a[0] << n) & (spint)0x7fffffffffffff;
}

// shift right by less than a word. Return shifted out part
static int modshr(unsigned int n, spint *a) {
  int i;
  spint r = a[0] & (((spint)1 << n) - (spint)1);
  for (i = 0; i < 6; i++) {
    a[i] = (a[i] >> n) | ((a[i + 1] << (55u - n)) & (spint)0x7fffffffffffff);
  }
  a[6] = a[6] >> n;
  return r;
}

// set a= 2^r
static void mod2r(unsigned int r, spint *a) {
  unsigned int n = r / 55u;
  unsigned int m = r % 55u;
  modzer(a);
  if (r >= 48 * 8)
    return;
  a[n] = 1;
  a[n] <<= m;
  nres(a, a);
}

// export to byte array
static void modexp(const spint *a, char *b) {
  int i;
  spint c[7];
  redc(a, c);
  for (i = 47; i >= 0; i--) {
    b[i] = c[0] & (spint)0xff;
    (void)modshr(8, c);
  }
}

// import from byte array
// returns 1 if in range, else 0
static int modimp(const char *b, spint *a) {
  int i, res;
  for (i = 0; i < 7; i++) {
    a[i] = 0;
  }
  for (i = 0; i < 48; i++) {
    modshl(8, a);
    a[0] += (spint)(unsigned char)b[i];
  }
  res = modfsb(a);
  nres(a, a);
  return res;
}

// determine sign
static int modsign(const spint *a) {
  spint c[7];
  redc(a, c);
  return c[0] % 2;
}

// return true if equal
static int modcmp(const spint *a, const spint *b) {
  spint c[7], d[7];
  int i, eq = 1;
  redc(a, c);
  redc(b, d);
  for (i = 0; i < 7; i++) {
    eq &= (((c[i] ^ d[i]) - 1) >> 55) & 1;
  }
  return eq;
}

// clang-format on
/******************************************************************************
 API functions calling generated code above
 ******************************************************************************/

#include <fp.h>

const digit_t ZERO[NWORDS_FIELD] = { 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0 };
const digit_t ONE[NWORDS_FIELD] = { 0x0000000000000007, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
                                    0x0000000000000000, 0x0000000000000000, 0x000e400000000000 };
// Montgomery representation of 2^-1
static const digit_t TWO_INV[NWORDS_FIELD] = { 0x0000000000000003, 0x0000000000000000, 0x0000000000000000,
                                               0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
                                               0x000f400000000000 };
// Montgomery representation of 3^-1
static const digit_t THREE_INV[NWORDS_FIELD] = { 0x0055555555555557, 0x002aaaaaaaaaaaaa, 0x0055555555555555,
                                                 0x002aaaaaaaaaaaaa, 0x0055555555555555, 0x002aaaaaaaaaaaaa,
                                                 0x000f955555555555 };
// Montgomery representation of 2^384
static const digit_t R2[NWORDS_FIELD] = { 0x0007e07e07e07e26, 0x007c0fc0fc0fc0fc, 0x0001f81f81f81f81,
                                          0x003f03f03f03f03f, 0x00607e07e07e07e0, 0x000fc0fc0fc0fc0f,
                                          0x000e9f81f81f81f8 };

void
fp_set_small(fp_t *x, const digit_t val)
{
    modint((int)val, *x);
}

void
fp_mul_small(fp_t *x, const fp_t *a, const uint32_t val)
{
    modmli(*a, (int)val, *x);
}

void
fp_set_zero(fp_t *x)
{
    modzer(*x);
}

void
fp_set_one(fp_t *x)
{
    modone(*x);
}

uint32_t
fp_is_equal(const fp_t *a, const fp_t *b)
{
    return -(uint32_t)modcmp(*a, *b);
}

uint32_t
fp_is_zero(const fp_t *a)
{
    return -(uint32_t)modis0(*a);
}

void
fp_copy(fp_t *out, const fp_t *a)
{
    modcpy(*a, *out);
}

void
fp_cswap(fp_t *a, fp_t *b, uint32_t ctl)
{
    modcsw((int)(ctl & 0x1), *a, *b);
}

void
fp_add(fp_t *out, const fp_t *a, const fp_t *b)
{
    modadd(*a, *b, *out);
}

void
fp_sub(fp_t *out, const fp_t *a, const fp_t *b)
{
    modsub(*a, *b, *out);
}

void
fp_neg(fp_t *out, const fp_t *a)
{
    modneg(*a, *out);
}

void
fp_sqr(fp_t *out, const fp_t *a)
{
    modsqr(*a, *out);
}

void
fp_mul(fp_t *out, const fp_t *a, const fp_t *b)
{
    modmul(*a, *b, *out);
}

void
fp_inv(fp_t *x)
{
    modinv(*x, NULL, *x);
}

uint32_t
fp_is_square(const fp_t *a)
{
    return -(uint32_t)modqr(NULL, *a);
}

void
fp_sqrt(fp_t *a)
{
    modsqrt(*a, NULL, *a);
}

void
fp_half(fp_t *out, const fp_t *a)
{
    modmul(TWO_INV, *a, *out);
}

void
fp_exp3div4(fp_t *out, const fp_t *a)
{
    modpro(*a, *out);
}

void
fp_div3(fp_t *out, const fp_t *a)
{
    modmul(THREE_INV, *a, *out);
}

void
fp_encode(void *dst, const fp_t *a)
{
    // Modified version of modexp()
    int i;
    spint c[7];
    redc(*a, c);
    for (i = 0; i < 48; i++) {
        ((char *)dst)[i] = c[0] & (spint)0xff;
        (void)modshr(8, c);
    }
}

uint32_t
fp_decode(fp_t *d, const void *src)
{
    // Modified version of modimp()
    int i;
    spint res;
    const unsigned char *b = src;
    for (i = 0; i < 7; i++) {
        (*d)[i] = 0;
    }
    for (i = 47; i >= 0; i--) {
        modshl(8, *d);
        (*d)[0] += (spint)b[i];
    }
    res = (spint)-modfsb(*d);
    nres(*d, *d);
    // If the value was canonical then res = -1; otherwise, res = 0
    for (i = 0; i < 7; i++) {
        (*d)[i] &= res;
    }
    return (uint32_t)res;
}

static inline unsigned char
add_carry(unsigned char cc, spint a, spint b, spint *d)
{
    udpint t = (udpint)a + (udpint)b + cc;
    *d = (spint)t;
    return (unsigned char)(t >> Wordlength);
}

static void
partial_reduce(spint *out, const spint *src)
{
    spint h, l, quo, rem;
    unsigned char cc;

    // Split value in high (8 bits) and low (376 bits) parts.
    h = src[5] >> 56;
    l = src[5] & 0x00FFFFFFFFFFFFFF;

    // 65*2^376 = 1 mod q; hence, we add floor(h/65) + (h mod 65)*2^376
    // to the low part.
    quo = (h * 0xFC1) >> 18;
    rem = h - (65 * quo);
    cc = add_carry(0, src[0], quo, &out[0]);
    cc = add_carry(cc, src[1], 0, &out[1]);
    cc = add_carry(cc, src[2], 0, &out[2]);
    cc = add_carry(cc, src[3], 0, &out[3]);
    cc = add_carry(cc, src[4], 0, &out[4]);
    (void)add_carry(cc, l, rem << 56, &out[5]);
}

// Little-endian encoding of a 64-bit integer.
static inline void
enc64le(void *dst, uint64_t x)
{
    uint8_t *buf = dst;
    buf[0] = (uint8_t)x;
    buf[1] = (uint8_t)(x >> 8);
    buf[2] = (uint8_t)(x >> 16);
    buf[3] = (uint8_t)(x >> 24);
    buf[4] = (uint8_t)(x >> 32);
    buf[5] = (uint8_t)(x >> 40);
    buf[6] = (uint8_t)(x >> 48);
    buf[7] = (uint8_t)(x >> 56);
}

// Little-endian decoding of a 64-bit integer.
static inline uint64_t
dec64le(const void *src)
{
    const uint8_t *buf = src;
    return (spint)buf[0] | ((spint)buf[1] << 8) | ((spint)buf[2] << 16) | ((spint)buf[3] << 24) |
           ((spint)buf[4] << 32) | ((spint)buf[5] << 40) | ((spint)buf[6] << 48) | ((spint)buf[7] << 56);
}

void
fp_decode_reduce(fp_t *d, const void *src, size_t len)
{
    uint64_t t[6];   // Stores Nbytes * 8 bits
    uint8_t tmp[48]; // Nbytes
    const uint8_t *b = src;

    fp_set_zero(d);
    if (len == 0) {
        return;
    }

    size_t rem = len % 48;
    if (rem != 0) {
        // Input size is not a multiple of 48, we decode a partial
        // block, which is already less than 2^376.
        size_t k = len - rem;
        memcpy(tmp, b + k, len - k);
        memset(tmp + len - k, 0, (sizeof tmp) - (len - k));
        fp_decode(d, tmp);
        len = k;
    }
    // Process all remaining blocks, in descending address order.
    while (len > 0) {
        fp_mul(d, d, &R2);
        len -= 48;
        t[0] = dec64le(b + len);
        t[1] = dec64le(b + len + 8);
        t[2] = dec64le(b + len + 16);
        t[3] = dec64le(b + len + 24);
        t[4] = dec64le(b + len + 32);
        t[5] = dec64le(b + len + 40);
        partial_reduce(t, t);
        enc64le(tmp, t[0]);
        enc64le(tmp + 8, t[1]);
        enc64le(tmp + 16, t[2]);
        enc64le(tmp + 24, t[3]);
        enc64le(tmp + 32, t[4]);
        enc64le(tmp + 40, t[5]);
        fp_t a;
        fp_decode(&a, tmp);
        fp_add(d, d, &a);
    }
}

// Vectorization
void prop_2(uint32x4_t *n) {
  __prop32(n);
}

// in: a, out: a/65
void divn(uint32x4_t* least, uint32x4_t* in){
  __div65(least, in);
}

void fp_bactched_reduction(uint32x4_t *out){
    uint32x4_t tmp;
    prop_2(out);
    divn(&tmp, out+(FP_LIMBS-1));
    out[0] = vaddq_u32(out[0], tmp);
}

void fp_add_batched(uint32x4_t* out, uint32x4_t *a, uint32x4_t *b){
    for (int i=0; i<FP_LIMBS; i++){
        out[i] = vaddq_u32(a[i], b[i]);
    }
}

void fp_sub_batched(uint32x4_t* out, uint32x4_t *a, uint32x4_t *b){
    uint32x4_t q[FP_LIMBS];
    for(int i = 0; i<(FP_LIMBS-1); i++) q[i] = vdupq_n_u32(Q2_VALUE);
    q[FP_LIMBS-1] = vdupq_n_u32(Q2_VALUE_SIGNIFICANT);

    for (int i=0; i<FP_LIMBS; i++){
        out[i] = vaddq_u32(a[i], vsubq_u32(q[i], b[i]));
    }
}

void fp_mul_batched(uint32x2_t *out, uint32x4_t *a, uint32x4_t *b){
    uint64x2_t t[2], mask = vdupq_n_u64(((uint64_t)1<<PER_LIMB)-1);
    uint32x2_t mod = vdup_n_u32(MASK_VALUE);
    uint32x2_t tmp[FP_LIMBS*2];
    
    //---x^0----
    t[0] = vmull_u32(((uint32x2_t*)a)[0], ((uint32x2_t*)b)[0]);
    t[1] = vmull_high_u32(a[0], b[0]);
    tmp[0] = vmovn_u64( vandq_u64(t[0], mask));
    tmp[1] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], 28);
    t[1] = vshrq_n_u64(t[1], 28);
    //---x^1----
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[0], ((uint32x2_t*)b)[2]);
    t[1] = vmlal_high_u32(t[1], a[0], b[1]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[2], ((uint32x2_t*)b)[0]);
    t[1] = vmlal_high_u32(t[1], a[1], b[0]);
    tmp[2] = vmovn_u64( vandq_u64(t[0], mask));
    tmp[3] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], 28);
    t[1] = vshrq_n_u64(t[1], 28);
    //---x^2----
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[0], ((uint32x2_t*)b)[4]);
    t[1] = vmlal_high_u32(t[1], a[0], b[2]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[2], ((uint32x2_t*)b)[2]);
    t[1] = vmlal_high_u32(t[1], a[1], b[1]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[4], ((uint32x2_t*)b)[0]);
    t[1] = vmlal_high_u32(t[1], a[2], b[0]);
    tmp[4] = vmovn_u64( vandq_u64(t[0], mask));
    tmp[5] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], 28);
    t[1] = vshrq_n_u64(t[1], 28);
    //---x^3----
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[0], ((uint32x2_t*)b)[6]);
    t[1] = vmlal_high_u32(t[1], a[0], b[3]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[2], ((uint32x2_t*)b)[4]);
    t[1] = vmlal_high_u32(t[1], a[1], b[2]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[4], ((uint32x2_t*)b)[2]);
    t[1] = vmlal_high_u32(t[1], a[2], b[1]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[6], ((uint32x2_t*)b)[0]);
    t[1] = vmlal_high_u32(t[1], a[3], b[0]);
    tmp[6] = vmovn_u64( vandq_u64(t[0], mask));
    tmp[7] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], 28);
    t[1] = vshrq_n_u64(t[1], 28);
    //---x^4----
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[0], ((uint32x2_t*)b)[8]);
    t[1] = vmlal_high_u32(t[1], a[0], b[4]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[2], ((uint32x2_t*)b)[6]);
    t[1] = vmlal_high_u32(t[1], a[1], b[3]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[4], ((uint32x2_t*)b)[4]);
    t[1] = vmlal_high_u32(t[1], a[2], b[2]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[6], ((uint32x2_t*)b)[2]);
    t[1] = vmlal_high_u32(t[1], a[3], b[1]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[8], ((uint32x2_t*)b)[0]);
    t[1] = vmlal_high_u32(t[1], a[4], b[0]);
    tmp[8] = vmovn_u64( vandq_u64(t[0], mask));
    tmp[9] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], 28);
    t[1] = vshrq_n_u64(t[1], 28);
    //---x^5----
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[0], ((uint32x2_t*)b)[10]);
    t[1] = vmlal_high_u32(t[1], a[0], b[5]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[2], ((uint32x2_t*)b)[8]);
    t[1] = vmlal_high_u32(t[1], a[1], b[4]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[4], ((uint32x2_t*)b)[6]);
    t[1] = vmlal_high_u32(t[1], a[2], b[3]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[6], ((uint32x2_t*)b)[4]);
    t[1] = vmlal_high_u32(t[1], a[3], b[2]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[8], ((uint32x2_t*)b)[2]);
    t[1] = vmlal_high_u32(t[1], a[4], b[1]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[10], ((uint32x2_t*)b)[0]);
    t[1] = vmlal_high_u32(t[1], a[5], b[0]);
    tmp[10] = vmovn_u64( vandq_u64(t[0], mask));
    tmp[11] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], 28);
    t[1] = vshrq_n_u64(t[1], 28);
    //---x^6----
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[0], ((uint32x2_t*)b)[12]);
    t[1] = vmlal_high_u32(t[1], a[0], b[6]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[2], ((uint32x2_t*)b)[10]);
    t[1] = vmlal_high_u32(t[1], a[1], b[5]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[4], ((uint32x2_t*)b)[8]);
    t[1] = vmlal_high_u32(t[1], a[2], b[4]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[6], ((uint32x2_t*)b)[6]);
    t[1] = vmlal_high_u32(t[1], a[3], b[3]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[8], ((uint32x2_t*)b)[4]);
    t[1] = vmlal_high_u32(t[1], a[4], b[2]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[10], ((uint32x2_t*)b)[2]);
    t[1] = vmlal_high_u32(t[1], a[5], b[1]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[12], ((uint32x2_t*)b)[0]);
    t[1] = vmlal_high_u32(t[1], a[6], b[0]);
    tmp[12] = vmovn_u64( vandq_u64(t[0], mask));
    tmp[13] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], 28);
    t[1] = vshrq_n_u64(t[1], 28);
    //---x^7----
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[0], ((uint32x2_t*)b)[14]);
    t[1] = vmlal_high_u32(t[1], a[0], b[7]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[2], ((uint32x2_t*)b)[12]);
    t[1] = vmlal_high_u32(t[1], a[1], b[6]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[4], ((uint32x2_t*)b)[10]);
    t[1] = vmlal_high_u32(t[1], a[2], b[5]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[6], ((uint32x2_t*)b)[8]);
    t[1] = vmlal_high_u32(t[1], a[3], b[4]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[8], ((uint32x2_t*)b)[6]);
    t[1] = vmlal_high_u32(t[1], a[4], b[3]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[10], ((uint32x2_t*)b)[4]);
    t[1] = vmlal_high_u32(t[1], a[5], b[2]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[12], ((uint32x2_t*)b)[2]);
    t[1] = vmlal_high_u32(t[1], a[6], b[1]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[14], ((uint32x2_t*)b)[0]);
    t[1] = vmlal_high_u32(t[1], a[7], b[0]);
    tmp[14] = vmovn_u64( vandq_u64(t[0], mask));
    tmp[15] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], 28);
    t[1] = vshrq_n_u64(t[1], 28);
    //---x^8----
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[0], ((uint32x2_t*)b)[16]);
    t[1] = vmlal_high_u32(t[1], a[0], b[8]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[2], ((uint32x2_t*)b)[14]);
    t[1] = vmlal_high_u32(t[1], a[1], b[7]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[4], ((uint32x2_t*)b)[12]);
    t[1] = vmlal_high_u32(t[1], a[2], b[6]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[6], ((uint32x2_t*)b)[10]);
    t[1] = vmlal_high_u32(t[1], a[3], b[5]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[8], ((uint32x2_t*)b)[8]);
    t[1] = vmlal_high_u32(t[1], a[4], b[4]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[10], ((uint32x2_t*)b)[6]);
    t[1] = vmlal_high_u32(t[1], a[5], b[3]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[12], ((uint32x2_t*)b)[4]);
    t[1] = vmlal_high_u32(t[1], a[6], b[2]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[14], ((uint32x2_t*)b)[2]);
    t[1] = vmlal_high_u32(t[1], a[7], b[1]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[16], ((uint32x2_t*)b)[0]);
    t[1] = vmlal_high_u32(t[1], a[8], b[0]);
    tmp[16] = vmovn_u64( vandq_u64(t[0], mask));
    tmp[17] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], 28);
    t[1] = vshrq_n_u64(t[1], 28);
    //---x^9----
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[0], ((uint32x2_t*)b)[18]);
    t[1] = vmlal_high_u32(t[1], a[0], b[9]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[2], ((uint32x2_t*)b)[16]);
    t[1] = vmlal_high_u32(t[1], a[1], b[8]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[4], ((uint32x2_t*)b)[14]);
    t[1] = vmlal_high_u32(t[1], a[2], b[7]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[6], ((uint32x2_t*)b)[12]);
    t[1] = vmlal_high_u32(t[1], a[3], b[6]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[8], ((uint32x2_t*)b)[10]);
    t[1] = vmlal_high_u32(t[1], a[4], b[5]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[10], ((uint32x2_t*)b)[8]);
    t[1] = vmlal_high_u32(t[1], a[5], b[4]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[12], ((uint32x2_t*)b)[6]);
    t[1] = vmlal_high_u32(t[1], a[6], b[3]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[14], ((uint32x2_t*)b)[4]);
    t[1] = vmlal_high_u32(t[1], a[7], b[2]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[16], ((uint32x2_t*)b)[2]);
    t[1] = vmlal_high_u32(t[1], a[8], b[1]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[18], ((uint32x2_t*)b)[0]);
    t[1] = vmlal_high_u32(t[1], a[9], b[0]);
    tmp[18] = vmovn_u64( vandq_u64(t[0], mask));
    tmp[19] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], 28);
    t[1] = vshrq_n_u64(t[1], 28);
    //---x^10----
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[0], ((uint32x2_t*)b)[20]);
    t[1] = vmlal_high_u32(t[1], a[0], b[10]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[2], ((uint32x2_t*)b)[18]);
    t[1] = vmlal_high_u32(t[1], a[1], b[9]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[4], ((uint32x2_t*)b)[16]);
    t[1] = vmlal_high_u32(t[1], a[2], b[8]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[6], ((uint32x2_t*)b)[14]);
    t[1] = vmlal_high_u32(t[1], a[3], b[7]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[8], ((uint32x2_t*)b)[12]);
    t[1] = vmlal_high_u32(t[1], a[4], b[6]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[10], ((uint32x2_t*)b)[10]);
    t[1] = vmlal_high_u32(t[1], a[5], b[5]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[12], ((uint32x2_t*)b)[8]);
    t[1] = vmlal_high_u32(t[1], a[6], b[4]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[14], ((uint32x2_t*)b)[6]);
    t[1] = vmlal_high_u32(t[1], a[7], b[3]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[16], ((uint32x2_t*)b)[4]);
    t[1] = vmlal_high_u32(t[1], a[8], b[2]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[18], ((uint32x2_t*)b)[2]);
    t[1] = vmlal_high_u32(t[1], a[9], b[1]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[20], ((uint32x2_t*)b)[0]);
    t[1] = vmlal_high_u32(t[1], a[10], b[0]);
    tmp[20] = vmovn_u64( vandq_u64(t[0], mask));
    tmp[21] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], 28);
    t[1] = vshrq_n_u64(t[1], 28);
    //---x^11----
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[0], ((uint32x2_t*)b)[22]);
    t[1] = vmlal_high_u32(t[1], a[0], b[11]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[2], ((uint32x2_t*)b)[20]);
    t[1] = vmlal_high_u32(t[1], a[1], b[10]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[4], ((uint32x2_t*)b)[18]);
    t[1] = vmlal_high_u32(t[1], a[2], b[9]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[6], ((uint32x2_t*)b)[16]);
    t[1] = vmlal_high_u32(t[1], a[3], b[8]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[8], ((uint32x2_t*)b)[14]);
    t[1] = vmlal_high_u32(t[1], a[4], b[7]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[10], ((uint32x2_t*)b)[12]);
    t[1] = vmlal_high_u32(t[1], a[5], b[6]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[12], ((uint32x2_t*)b)[10]);
    t[1] = vmlal_high_u32(t[1], a[6], b[5]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[14], ((uint32x2_t*)b)[8]);
    t[1] = vmlal_high_u32(t[1], a[7], b[4]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[16], ((uint32x2_t*)b)[6]);
    t[1] = vmlal_high_u32(t[1], a[8], b[3]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[18], ((uint32x2_t*)b)[4]);
    t[1] = vmlal_high_u32(t[1], a[9], b[2]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[20], ((uint32x2_t*)b)[2]);
    t[1] = vmlal_high_u32(t[1], a[10], b[1]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[22], ((uint32x2_t*)b)[0]);
    t[1] = vmlal_high_u32(t[1], a[11], b[0]);
    tmp[22] = vmovn_u64( vandq_u64(t[0], mask));
    tmp[23] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], 28);
    t[1] = vshrq_n_u64(t[1], 28);
    //---x^12----
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[0], ((uint32x2_t*)b)[24]);
    t[1] = vmlal_high_u32(t[1], a[0], b[12]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[2], ((uint32x2_t*)b)[22]);
    t[1] = vmlal_high_u32(t[1], a[1], b[11]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[4], ((uint32x2_t*)b)[20]);
    t[1] = vmlal_high_u32(t[1], a[2], b[10]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[6], ((uint32x2_t*)b)[18]);
    t[1] = vmlal_high_u32(t[1], a[3], b[9]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[8], ((uint32x2_t*)b)[16]);
    t[1] = vmlal_high_u32(t[1], a[4], b[8]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[10], ((uint32x2_t*)b)[14]);
    t[1] = vmlal_high_u32(t[1], a[5], b[7]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[12], ((uint32x2_t*)b)[12]);
    t[1] = vmlal_high_u32(t[1], a[6], b[6]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[14], ((uint32x2_t*)b)[10]);
    t[1] = vmlal_high_u32(t[1], a[7], b[5]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[16], ((uint32x2_t*)b)[8]);
    t[1] = vmlal_high_u32(t[1], a[8], b[4]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[18], ((uint32x2_t*)b)[6]);
    t[1] = vmlal_high_u32(t[1], a[9], b[3]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[20], ((uint32x2_t*)b)[4]);
    t[1] = vmlal_high_u32(t[1], a[10], b[2]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[22], ((uint32x2_t*)b)[2]);
    t[1] = vmlal_high_u32(t[1], a[11], b[1]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[24], ((uint32x2_t*)b)[0]);
    t[1] = vmlal_high_u32(t[1], a[12], b[0]);
    tmp[24] = vmovn_u64( vandq_u64(t[0], mask));
    tmp[25] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], 28);
    t[1] = vshrq_n_u64(t[1], 28);
    //---x^13----
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[0], ((uint32x2_t*)b)[26]);
    t[1] = vmlal_high_u32(t[1], a[0], b[13]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[2], ((uint32x2_t*)b)[24]);
    t[1] = vmlal_high_u32(t[1], a[1], b[12]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[4], ((uint32x2_t*)b)[22]);
    t[1] = vmlal_high_u32(t[1], a[2], b[11]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[6], ((uint32x2_t*)b)[20]);
    t[1] = vmlal_high_u32(t[1], a[3], b[10]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[8], ((uint32x2_t*)b)[18]);
    t[1] = vmlal_high_u32(t[1], a[4], b[9]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[10], ((uint32x2_t*)b)[16]);
    t[1] = vmlal_high_u32(t[1], a[5], b[8]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[12], ((uint32x2_t*)b)[14]);
    t[1] = vmlal_high_u32(t[1], a[6], b[7]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[14], ((uint32x2_t*)b)[12]);
    t[1] = vmlal_high_u32(t[1], a[7], b[6]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[16], ((uint32x2_t*)b)[10]);
    t[1] = vmlal_high_u32(t[1], a[8], b[5]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[18], ((uint32x2_t*)b)[8]);
    t[1] = vmlal_high_u32(t[1], a[9], b[4]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[20], ((uint32x2_t*)b)[6]);
    t[1] = vmlal_high_u32(t[1], a[10], b[3]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[22], ((uint32x2_t*)b)[4]);
    t[1] = vmlal_high_u32(t[1], a[11], b[2]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[24], ((uint32x2_t*)b)[2]);
    t[1] = vmlal_high_u32(t[1], a[12], b[1]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[26], ((uint32x2_t*)b)[0]);
    t[1] = vmlal_high_u32(t[1], a[13], b[0]);
    t[0] = vmlal_u32(t[0], tmp[0], mod);
    t[1] = vmlal_u32(t[1], tmp[1], mod);
    tmp[26] = vmovn_u64( vandq_u64(t[0], mask));
    tmp[27] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], 28);
    t[1] = vshrq_n_u64(t[1], 28);
    //---x^14----
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[2], ((uint32x2_t*)b)[26]);
    t[1] = vmlal_high_u32(t[1], a[1], b[13]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[4], ((uint32x2_t*)b)[24]);
    t[1] = vmlal_high_u32(t[1], a[2], b[12]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[6], ((uint32x2_t*)b)[22]);
    t[1] = vmlal_high_u32(t[1], a[3], b[11]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[8], ((uint32x2_t*)b)[20]);
    t[1] = vmlal_high_u32(t[1], a[4], b[10]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[10], ((uint32x2_t*)b)[18]);
    t[1] = vmlal_high_u32(t[1], a[5], b[9]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[12], ((uint32x2_t*)b)[16]);
    t[1] = vmlal_high_u32(t[1], a[6], b[8]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[14], ((uint32x2_t*)b)[14]);
    t[1] = vmlal_high_u32(t[1], a[7], b[7]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[16], ((uint32x2_t*)b)[12]);
    t[1] = vmlal_high_u32(t[1], a[8], b[6]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[18], ((uint32x2_t*)b)[10]);
    t[1] = vmlal_high_u32(t[1], a[9], b[5]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[20], ((uint32x2_t*)b)[8]);
    t[1] = vmlal_high_u32(t[1], a[10], b[4]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[22], ((uint32x2_t*)b)[6]);
    t[1] = vmlal_high_u32(t[1], a[11], b[3]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[24], ((uint32x2_t*)b)[4]);
    t[1] = vmlal_high_u32(t[1], a[12], b[2]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[26], ((uint32x2_t*)b)[2]);
    t[1] = vmlal_high_u32(t[1], a[13], b[1]);
    t[0] = vmlal_u32(t[0], tmp[2], mod);
    t[1] = vmlal_u32(t[1], tmp[3], mod);
    out[0] = vmovn_u64( vandq_u64(t[0], mask));
    out[1] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], 28);
    t[1] = vshrq_n_u64(t[1], 28);
    //---x^15----
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[4], ((uint32x2_t*)b)[26]);
    t[1] = vmlal_high_u32(t[1], a[2], b[13]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[6], ((uint32x2_t*)b)[24]);
    t[1] = vmlal_high_u32(t[1], a[3], b[12]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[8], ((uint32x2_t*)b)[22]);
    t[1] = vmlal_high_u32(t[1], a[4], b[11]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[10], ((uint32x2_t*)b)[20]);
    t[1] = vmlal_high_u32(t[1], a[5], b[10]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[12], ((uint32x2_t*)b)[18]);
    t[1] = vmlal_high_u32(t[1], a[6], b[9]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[14], ((uint32x2_t*)b)[16]);
    t[1] = vmlal_high_u32(t[1], a[7], b[8]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[16], ((uint32x2_t*)b)[14]);
    t[1] = vmlal_high_u32(t[1], a[8], b[7]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[18], ((uint32x2_t*)b)[12]);
    t[1] = vmlal_high_u32(t[1], a[9], b[6]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[20], ((uint32x2_t*)b)[10]);
    t[1] = vmlal_high_u32(t[1], a[10], b[5]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[22], ((uint32x2_t*)b)[8]);
    t[1] = vmlal_high_u32(t[1], a[11], b[4]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[24], ((uint32x2_t*)b)[6]);
    t[1] = vmlal_high_u32(t[1], a[12], b[3]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[26], ((uint32x2_t*)b)[4]);
    t[1] = vmlal_high_u32(t[1], a[13], b[2]);
    t[0] = vmlal_u32(t[0], tmp[4], mod);
    t[1] = vmlal_u32(t[1], tmp[5], mod);
    out[2] = vmovn_u64( vandq_u64(t[0], mask));
    out[3] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], 28);
    t[1] = vshrq_n_u64(t[1], 28);
    //---x^16----
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[6], ((uint32x2_t*)b)[26]);
    t[1] = vmlal_high_u32(t[1], a[3], b[13]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[8], ((uint32x2_t*)b)[24]);
    t[1] = vmlal_high_u32(t[1], a[4], b[12]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[10], ((uint32x2_t*)b)[22]);
    t[1] = vmlal_high_u32(t[1], a[5], b[11]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[12], ((uint32x2_t*)b)[20]);
    t[1] = vmlal_high_u32(t[1], a[6], b[10]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[14], ((uint32x2_t*)b)[18]);
    t[1] = vmlal_high_u32(t[1], a[7], b[9]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[16], ((uint32x2_t*)b)[16]);
    t[1] = vmlal_high_u32(t[1], a[8], b[8]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[18], ((uint32x2_t*)b)[14]);
    t[1] = vmlal_high_u32(t[1], a[9], b[7]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[20], ((uint32x2_t*)b)[12]);
    t[1] = vmlal_high_u32(t[1], a[10], b[6]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[22], ((uint32x2_t*)b)[10]);
    t[1] = vmlal_high_u32(t[1], a[11], b[5]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[24], ((uint32x2_t*)b)[8]);
    t[1] = vmlal_high_u32(t[1], a[12], b[4]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[26], ((uint32x2_t*)b)[6]);
    t[1] = vmlal_high_u32(t[1], a[13], b[3]);
    t[0] = vmlal_u32(t[0], tmp[6], mod);
    t[1] = vmlal_u32(t[1], tmp[7], mod);
    out[4] = vmovn_u64( vandq_u64(t[0], mask));
    out[5] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], 28);
    t[1] = vshrq_n_u64(t[1], 28);
    //---x^17----
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[8], ((uint32x2_t*)b)[26]);
    t[1] = vmlal_high_u32(t[1], a[4], b[13]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[10], ((uint32x2_t*)b)[24]);
    t[1] = vmlal_high_u32(t[1], a[5], b[12]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[12], ((uint32x2_t*)b)[22]);
    t[1] = vmlal_high_u32(t[1], a[6], b[11]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[14], ((uint32x2_t*)b)[20]);
    t[1] = vmlal_high_u32(t[1], a[7], b[10]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[16], ((uint32x2_t*)b)[18]);
    t[1] = vmlal_high_u32(t[1], a[8], b[9]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[18], ((uint32x2_t*)b)[16]);
    t[1] = vmlal_high_u32(t[1], a[9], b[8]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[20], ((uint32x2_t*)b)[14]);
    t[1] = vmlal_high_u32(t[1], a[10], b[7]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[22], ((uint32x2_t*)b)[12]);
    t[1] = vmlal_high_u32(t[1], a[11], b[6]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[24], ((uint32x2_t*)b)[10]);
    t[1] = vmlal_high_u32(t[1], a[12], b[5]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[26], ((uint32x2_t*)b)[8]);
    t[1] = vmlal_high_u32(t[1], a[13], b[4]);
    t[0] = vmlal_u32(t[0], tmp[8], mod);
    t[1] = vmlal_u32(t[1], tmp[9], mod);
    out[6] = vmovn_u64( vandq_u64(t[0], mask));
    out[7] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], 28);
    t[1] = vshrq_n_u64(t[1], 28);
    //---x^18----
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[10], ((uint32x2_t*)b)[26]);
    t[1] = vmlal_high_u32(t[1], a[5], b[13]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[12], ((uint32x2_t*)b)[24]);
    t[1] = vmlal_high_u32(t[1], a[6], b[12]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[14], ((uint32x2_t*)b)[22]);
    t[1] = vmlal_high_u32(t[1], a[7], b[11]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[16], ((uint32x2_t*)b)[20]);
    t[1] = vmlal_high_u32(t[1], a[8], b[10]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[18], ((uint32x2_t*)b)[18]);
    t[1] = vmlal_high_u32(t[1], a[9], b[9]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[20], ((uint32x2_t*)b)[16]);
    t[1] = vmlal_high_u32(t[1], a[10], b[8]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[22], ((uint32x2_t*)b)[14]);
    t[1] = vmlal_high_u32(t[1], a[11], b[7]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[24], ((uint32x2_t*)b)[12]);
    t[1] = vmlal_high_u32(t[1], a[12], b[6]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[26], ((uint32x2_t*)b)[10]);
    t[1] = vmlal_high_u32(t[1], a[13], b[5]);
    t[0] = vmlal_u32(t[0], tmp[10], mod);
    t[1] = vmlal_u32(t[1], tmp[11], mod);
    out[8] = vmovn_u64( vandq_u64(t[0], mask));
    out[9] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], 28);
    t[1] = vshrq_n_u64(t[1], 28);
    //---x^19----
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[12], ((uint32x2_t*)b)[26]);
    t[1] = vmlal_high_u32(t[1], a[6], b[13]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[14], ((uint32x2_t*)b)[24]);
    t[1] = vmlal_high_u32(t[1], a[7], b[12]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[16], ((uint32x2_t*)b)[22]);
    t[1] = vmlal_high_u32(t[1], a[8], b[11]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[18], ((uint32x2_t*)b)[20]);
    t[1] = vmlal_high_u32(t[1], a[9], b[10]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[20], ((uint32x2_t*)b)[18]);
    t[1] = vmlal_high_u32(t[1], a[10], b[9]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[22], ((uint32x2_t*)b)[16]);
    t[1] = vmlal_high_u32(t[1], a[11], b[8]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[24], ((uint32x2_t*)b)[14]);
    t[1] = vmlal_high_u32(t[1], a[12], b[7]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[26], ((uint32x2_t*)b)[12]);
    t[1] = vmlal_high_u32(t[1], a[13], b[6]);
    t[0] = vmlal_u32(t[0], tmp[12], mod);
    t[1] = vmlal_u32(t[1], tmp[13], mod);
    out[10] = vmovn_u64( vandq_u64(t[0], mask));
    out[11] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], 28);
    t[1] = vshrq_n_u64(t[1], 28);
    //---x^20----
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[14], ((uint32x2_t*)b)[26]);
    t[1] = vmlal_high_u32(t[1], a[7], b[13]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[16], ((uint32x2_t*)b)[24]);
    t[1] = vmlal_high_u32(t[1], a[8], b[12]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[18], ((uint32x2_t*)b)[22]);
    t[1] = vmlal_high_u32(t[1], a[9], b[11]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[20], ((uint32x2_t*)b)[20]);
    t[1] = vmlal_high_u32(t[1], a[10], b[10]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[22], ((uint32x2_t*)b)[18]);
    t[1] = vmlal_high_u32(t[1], a[11], b[9]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[24], ((uint32x2_t*)b)[16]);
    t[1] = vmlal_high_u32(t[1], a[12], b[8]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[26], ((uint32x2_t*)b)[14]);
    t[1] = vmlal_high_u32(t[1], a[13], b[7]);
    t[0] = vmlal_u32(t[0], tmp[14], mod);
    t[1] = vmlal_u32(t[1], tmp[15], mod);
    out[12] = vmovn_u64( vandq_u64(t[0], mask));
    out[13] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], 28);
    t[1] = vshrq_n_u64(t[1], 28);
    //---x^21----
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[16], ((uint32x2_t*)b)[26]);
    t[1] = vmlal_high_u32(t[1], a[8], b[13]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[18], ((uint32x2_t*)b)[24]);
    t[1] = vmlal_high_u32(t[1], a[9], b[12]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[20], ((uint32x2_t*)b)[22]);
    t[1] = vmlal_high_u32(t[1], a[10], b[11]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[22], ((uint32x2_t*)b)[20]);
    t[1] = vmlal_high_u32(t[1], a[11], b[10]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[24], ((uint32x2_t*)b)[18]);
    t[1] = vmlal_high_u32(t[1], a[12], b[9]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[26], ((uint32x2_t*)b)[16]);
    t[1] = vmlal_high_u32(t[1], a[13], b[8]);
    t[0] = vmlal_u32(t[0], tmp[16], mod);
    t[1] = vmlal_u32(t[1], tmp[17], mod);
    out[14] = vmovn_u64( vandq_u64(t[0], mask));
    out[15] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], 28);
    t[1] = vshrq_n_u64(t[1], 28);
    //---x^22----
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[18], ((uint32x2_t*)b)[26]);
    t[1] = vmlal_high_u32(t[1], a[9], b[13]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[20], ((uint32x2_t*)b)[24]);
    t[1] = vmlal_high_u32(t[1], a[10], b[12]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[22], ((uint32x2_t*)b)[22]);
    t[1] = vmlal_high_u32(t[1], a[11], b[11]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[24], ((uint32x2_t*)b)[20]);
    t[1] = vmlal_high_u32(t[1], a[12], b[10]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[26], ((uint32x2_t*)b)[18]);
    t[1] = vmlal_high_u32(t[1], a[13], b[9]);
    t[0] = vmlal_u32(t[0], tmp[18], mod);
    t[1] = vmlal_u32(t[1], tmp[19], mod);
    out[16] = vmovn_u64( vandq_u64(t[0], mask));
    out[17] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], 28);
    t[1] = vshrq_n_u64(t[1], 28);
    //---x^23----
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[20], ((uint32x2_t*)b)[26]);
    t[1] = vmlal_high_u32(t[1], a[10], b[13]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[22], ((uint32x2_t*)b)[24]);
    t[1] = vmlal_high_u32(t[1], a[11], b[12]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[24], ((uint32x2_t*)b)[22]);
    t[1] = vmlal_high_u32(t[1], a[12], b[11]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[26], ((uint32x2_t*)b)[20]);
    t[1] = vmlal_high_u32(t[1], a[13], b[10]);
    t[0] = vmlal_u32(t[0], tmp[20], mod);
    t[1] = vmlal_u32(t[1], tmp[21], mod);
    out[18] = vmovn_u64( vandq_u64(t[0], mask));
    out[19] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], 28);
    t[1] = vshrq_n_u64(t[1], 28);
    //---x^24----
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[22], ((uint32x2_t*)b)[26]);
    t[1] = vmlal_high_u32(t[1], a[11], b[13]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[24], ((uint32x2_t*)b)[24]);
    t[1] = vmlal_high_u32(t[1], a[12], b[12]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[26], ((uint32x2_t*)b)[22]);
    t[1] = vmlal_high_u32(t[1], a[13], b[11]);
    t[0] = vmlal_u32(t[0], tmp[22], mod);
    t[1] = vmlal_u32(t[1], tmp[23], mod);
    out[20] = vmovn_u64( vandq_u64(t[0], mask));
    out[21] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], 28);
    t[1] = vshrq_n_u64(t[1], 28);
    //---x^25----
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[24], ((uint32x2_t*)b)[26]);
    t[1] = vmlal_high_u32(t[1], a[12], b[13]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[26], ((uint32x2_t*)b)[24]);
    t[1] = vmlal_high_u32(t[1], a[13], b[12]);
    t[0] = vmlal_u32(t[0], tmp[24], mod);
    t[1] = vmlal_u32(t[1], tmp[25], mod);
    out[22] = vmovn_u64( vandq_u64(t[0], mask));
    out[23] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], 28);
    t[1] = vshrq_n_u64(t[1], 28);
    //---x^26----
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[26], ((uint32x2_t*)b)[26]);
    t[1] = vmlal_high_u32(t[1], a[13], b[13]);
    t[0] = vmlal_u32(t[0], tmp[26], mod);
    t[1] = vmlal_u32(t[1], tmp[27], mod);
    out[24] = vmovn_u64( vandq_u64(t[0], mask));
    out[25] = vmovn_u64( vandq_u64(t[1], mask));
    //---x^27----
    out[26] = vshrn_n_u64(t[0], 28);
    out[27] = vshrn_n_u64(t[1], 28);
}

void fp_sqr_batched(uint32x2_t *out, uint32x4_t *a){
    uint64x2_t t[2], mask = vdupq_n_u64(((uint64_t)1<<PER_LIMB)-1);
    uint32x2_t mod = vdup_n_u32(MASK_VALUE);
    uint32x2_t tmp[FP_LIMBS*2];
    
    //---x^0----
    t[0] = vmull_u32(((uint32x2_t*)a)[0], ((uint32x2_t*)a)[0]);
    t[1] = vmull_high_u32(a[0], a[0]);
    tmp[0] = vmovn_u64( vandq_u64(t[0], mask));
    tmp[1] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], 28);
    t[1] = vshrq_n_u64(t[1], 28);
    //---x^1----
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[0], ((int32x2_t*)a)[2]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[0], ((int32x4_t*)a)[1]);
    tmp[2] = vmovn_u64( vandq_u64(t[0], mask));
    tmp[3] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], PER_LIMB);
    t[1] = vshrq_n_u64(t[1], PER_LIMB);
    //---x^2----
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[0], ((int32x2_t*)a)[4]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[0], ((int32x4_t*)a)[2]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[2], ((uint32x2_t*)a)[2]);
    t[1] = vmlal_high_u32(t[1], a[1], a[1]);
    tmp[4] = vmovn_u64( vandq_u64(t[0], mask));
    tmp[5] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], PER_LIMB);
    t[1] = vshrq_n_u64(t[1], PER_LIMB);
    //---x^3----
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[0], ((int32x2_t*)a)[6]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[0], ((int32x4_t*)a)[3]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[2], ((int32x2_t*)a)[4]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[1], ((int32x4_t*)a)[2]);
    tmp[6] = vmovn_u64( vandq_u64(t[0], mask));
    tmp[7] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], PER_LIMB);
    t[1] = vshrq_n_u64(t[1], PER_LIMB);
    //---x^4----
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[0], ((int32x2_t*)a)[8]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[0], ((int32x4_t*)a)[4]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[2], ((int32x2_t*)a)[6]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[1], ((int32x4_t*)a)[3]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[4], ((uint32x2_t*)a)[4]);
    t[1] = vmlal_high_u32(t[1], a[2], a[2]);
    tmp[8] = vmovn_u64( vandq_u64(t[0], mask));
    tmp[9] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], PER_LIMB);
    t[1] = vshrq_n_u64(t[1], PER_LIMB);
    //---x^5----
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[0], ((int32x2_t*)a)[10]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[0], ((int32x4_t*)a)[5]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[2], ((int32x2_t*)a)[8]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[1], ((int32x4_t*)a)[4]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[4], ((int32x2_t*)a)[6]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[2], ((int32x4_t*)a)[3]);
    tmp[10] = vmovn_u64( vandq_u64(t[0], mask));
    tmp[11] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], PER_LIMB);
    t[1] = vshrq_n_u64(t[1], PER_LIMB);
    //---x^6----
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[0], ((int32x2_t*)a)[12]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[0], ((int32x4_t*)a)[6]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[2], ((int32x2_t*)a)[10]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[1], ((int32x4_t*)a)[5]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[4], ((int32x2_t*)a)[8]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[2], ((int32x4_t*)a)[4]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[6], ((uint32x2_t*)a)[6]);
    t[1] = vmlal_high_u32(t[1], a[3], a[3]);
    tmp[12] = vmovn_u64( vandq_u64(t[0], mask));
    tmp[13] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], PER_LIMB);
    t[1] = vshrq_n_u64(t[1], PER_LIMB);
    //---x^7----
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[0], ((int32x2_t*)a)[14]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[0], ((int32x4_t*)a)[7]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[2], ((int32x2_t*)a)[12]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[1], ((int32x4_t*)a)[6]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[4], ((int32x2_t*)a)[10]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[2], ((int32x4_t*)a)[5]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[6], ((int32x2_t*)a)[8]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[3], ((int32x4_t*)a)[4]);
    tmp[14] = vmovn_u64( vandq_u64(t[0], mask));
    tmp[15] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], PER_LIMB);
    t[1] = vshrq_n_u64(t[1], PER_LIMB);
    //---x^8----
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[0], ((int32x2_t*)a)[16]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[0], ((int32x4_t*)a)[8]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[2], ((int32x2_t*)a)[14]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[1], ((int32x4_t*)a)[7]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[4], ((int32x2_t*)a)[12]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[2], ((int32x4_t*)a)[6]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[6], ((int32x2_t*)a)[10]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[3], ((int32x4_t*)a)[5]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[8], ((uint32x2_t*)a)[8]);
    t[1] = vmlal_high_u32(t[1], a[4], a[4]);
    tmp[16] = vmovn_u64( vandq_u64(t[0], mask));
    tmp[17] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], PER_LIMB);
    t[1] = vshrq_n_u64(t[1], PER_LIMB);
    //---x^9----
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[0], ((int32x2_t*)a)[18]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[0], ((int32x4_t*)a)[9]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[2], ((int32x2_t*)a)[16]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[1], ((int32x4_t*)a)[8]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[4], ((int32x2_t*)a)[14]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[2], ((int32x4_t*)a)[7]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[6], ((int32x2_t*)a)[12]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[3], ((int32x4_t*)a)[6]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[8], ((int32x2_t*)a)[10]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[4], ((int32x4_t*)a)[5]);
    tmp[18] = vmovn_u64( vandq_u64(t[0], mask));
    tmp[19] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], PER_LIMB);
    t[1] = vshrq_n_u64(t[1], PER_LIMB);
    //---x^10----
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[0], ((int32x2_t*)a)[20]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[0], ((int32x4_t*)a)[10]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[2], ((int32x2_t*)a)[18]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[1], ((int32x4_t*)a)[9]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[4], ((int32x2_t*)a)[16]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[2], ((int32x4_t*)a)[8]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[6], ((int32x2_t*)a)[14]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[3], ((int32x4_t*)a)[7]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[8], ((int32x2_t*)a)[12]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[4], ((int32x4_t*)a)[6]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[10], ((uint32x2_t*)a)[10]);
    t[1] = vmlal_high_u32(t[1], a[5], a[5]);
    tmp[20] = vmovn_u64( vandq_u64(t[0], mask));
    tmp[21] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], PER_LIMB);
    t[1] = vshrq_n_u64(t[1], PER_LIMB);
    //---x^11----
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[0], ((int32x2_t*)a)[22]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[0], ((int32x4_t*)a)[11]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[2], ((int32x2_t*)a)[20]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[1], ((int32x4_t*)a)[10]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[4], ((int32x2_t*)a)[18]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[2], ((int32x4_t*)a)[9]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[6], ((int32x2_t*)a)[16]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[3], ((int32x4_t*)a)[8]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[8], ((int32x2_t*)a)[14]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[4], ((int32x4_t*)a)[7]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[10], ((int32x2_t*)a)[12]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[5], ((int32x4_t*)a)[6]);
    tmp[22] = vmovn_u64( vandq_u64(t[0], mask));
    tmp[23] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], PER_LIMB);
    t[1] = vshrq_n_u64(t[1], PER_LIMB);
    //---x^12----
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[0], ((int32x2_t*)a)[24]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[0], ((int32x4_t*)a)[12]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[2], ((int32x2_t*)a)[22]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[1], ((int32x4_t*)a)[11]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[4], ((int32x2_t*)a)[20]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[2], ((int32x4_t*)a)[10]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[6], ((int32x2_t*)a)[18]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[3], ((int32x4_t*)a)[9]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[8], ((int32x2_t*)a)[16]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[4], ((int32x4_t*)a)[8]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[10], ((int32x2_t*)a)[14]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[5], ((int32x4_t*)a)[7]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[12], ((uint32x2_t*)a)[12]);
    t[1] = vmlal_high_u32(t[1], a[6], a[6]);
    tmp[24] = vmovn_u64( vandq_u64(t[0], mask));
    tmp[25] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], PER_LIMB);
    t[1] = vshrq_n_u64(t[1], PER_LIMB);
    //---x^13----
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[0], ((int32x2_t*)a)[26]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[0], ((int32x4_t*)a)[13]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[2], ((int32x2_t*)a)[24]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[1], ((int32x4_t*)a)[12]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[4], ((int32x2_t*)a)[22]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[2], ((int32x4_t*)a)[11]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[6], ((int32x2_t*)a)[20]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[3], ((int32x4_t*)a)[10]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[8], ((int32x2_t*)a)[18]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[4], ((int32x4_t*)a)[9]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[10], ((int32x2_t*)a)[16]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[5], ((int32x4_t*)a)[8]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[12], ((int32x2_t*)a)[14]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[6], ((int32x4_t*)a)[7]);
    t[0] = vmlal_u32(t[0], tmp[0], mod);
    t[1] = vmlal_u32(t[1], tmp[1], mod);
    tmp[26] = vmovn_u64( vandq_u64(t[0], mask));
    tmp[27] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], PER_LIMB);
    t[1] = vshrq_n_u64(t[1], PER_LIMB);
    //---x^14----
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[2], ((int32x2_t*)a)[26]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[1], ((int32x4_t*)a)[13]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[4], ((int32x2_t*)a)[24]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[2], ((int32x4_t*)a)[12]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[6], ((int32x2_t*)a)[22]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[3], ((int32x4_t*)a)[11]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[8], ((int32x2_t*)a)[20]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[4], ((int32x4_t*)a)[10]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[10], ((int32x2_t*)a)[18]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[5], ((int32x4_t*)a)[9]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[12], ((int32x2_t*)a)[16]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[6], ((int32x4_t*)a)[8]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[14], ((uint32x2_t*)a)[14]);
    t[1] = vmlal_high_u32(t[1], a[7], a[7]);
    t[0] = vmlal_u32(t[0], tmp[2], mod);
    t[1] = vmlal_u32(t[1], tmp[3], mod);
    out[0] = vmovn_u64( vandq_u64(t[0], mask));
    out[1] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], PER_LIMB);
    t[1] = vshrq_n_u64(t[1], PER_LIMB);
    //---x^15----
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[4], ((int32x2_t*)a)[26]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[2], ((int32x4_t*)a)[13]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[6], ((int32x2_t*)a)[24]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[3], ((int32x4_t*)a)[12]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[8], ((int32x2_t*)a)[22]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[4], ((int32x4_t*)a)[11]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[10], ((int32x2_t*)a)[20]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[5], ((int32x4_t*)a)[10]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[12], ((int32x2_t*)a)[18]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[6], ((int32x4_t*)a)[9]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[14], ((int32x2_t*)a)[16]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[7], ((int32x4_t*)a)[8]);
    t[0] = vmlal_u32(t[0], tmp[4], mod);
    t[1] = vmlal_u32(t[1], tmp[5], mod);
    out[2] = vmovn_u64( vandq_u64(t[0], mask));
    out[3] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], PER_LIMB);
    t[1] = vshrq_n_u64(t[1], PER_LIMB);
    //---x^16----
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[6], ((int32x2_t*)a)[26]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[3], ((int32x4_t*)a)[13]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[8], ((int32x2_t*)a)[24]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[4], ((int32x4_t*)a)[12]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[10], ((int32x2_t*)a)[22]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[5], ((int32x4_t*)a)[11]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[12], ((int32x2_t*)a)[20]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[6], ((int32x4_t*)a)[10]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[14], ((int32x2_t*)a)[18]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[7], ((int32x4_t*)a)[9]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[16], ((uint32x2_t*)a)[16]);
    t[1] = vmlal_high_u32(t[1], a[8], a[8]);
    t[0] = vmlal_u32(t[0], tmp[6], mod);
    t[1] = vmlal_u32(t[1], tmp[7], mod);
    out[4] = vmovn_u64( vandq_u64(t[0], mask));
    out[5] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], PER_LIMB);
    t[1] = vshrq_n_u64(t[1], PER_LIMB);
    //---x^17----
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[8], ((int32x2_t*)a)[26]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[4], ((int32x4_t*)a)[13]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[10], ((int32x2_t*)a)[24]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[5], ((int32x4_t*)a)[12]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[12], ((int32x2_t*)a)[22]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[6], ((int32x4_t*)a)[11]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[14], ((int32x2_t*)a)[20]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[7], ((int32x4_t*)a)[10]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[16], ((int32x2_t*)a)[18]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[8], ((int32x4_t*)a)[9]);
    t[0] = vmlal_u32(t[0], tmp[8], mod);
    t[1] = vmlal_u32(t[1], tmp[9], mod);
    out[6] = vmovn_u64( vandq_u64(t[0], mask));
    out[7] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], PER_LIMB);
    t[1] = vshrq_n_u64(t[1], PER_LIMB);
    //---x^18----
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[10], ((int32x2_t*)a)[26]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[5], ((int32x4_t*)a)[13]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[12], ((int32x2_t*)a)[24]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[6], ((int32x4_t*)a)[12]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[14], ((int32x2_t*)a)[22]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[7], ((int32x4_t*)a)[11]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[16], ((int32x2_t*)a)[20]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[8], ((int32x4_t*)a)[10]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[18], ((uint32x2_t*)a)[18]);
    t[1] = vmlal_high_u32(t[1], a[9], a[9]);
    t[0] = vmlal_u32(t[0], tmp[10], mod);
    t[1] = vmlal_u32(t[1], tmp[11], mod);
    out[8] = vmovn_u64( vandq_u64(t[0], mask));
    out[9] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], PER_LIMB);
    t[1] = vshrq_n_u64(t[1], PER_LIMB);
    //---x^19----
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[12], ((int32x2_t*)a)[26]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[6], ((int32x4_t*)a)[13]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[14], ((int32x2_t*)a)[24]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[7], ((int32x4_t*)a)[12]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[16], ((int32x2_t*)a)[22]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[8], ((int32x4_t*)a)[11]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[18], ((int32x2_t*)a)[20]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[9], ((int32x4_t*)a)[10]);
    t[0] = vmlal_u32(t[0], tmp[12], mod);
    t[1] = vmlal_u32(t[1], tmp[13], mod);
    out[10] = vmovn_u64( vandq_u64(t[0], mask));
    out[11] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], PER_LIMB);
    t[1] = vshrq_n_u64(t[1], PER_LIMB);
    //---x^20----
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[14], ((int32x2_t*)a)[26]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[7], ((int32x4_t*)a)[13]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[16], ((int32x2_t*)a)[24]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[8], ((int32x4_t*)a)[12]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[18], ((int32x2_t*)a)[22]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[9], ((int32x4_t*)a)[11]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[20], ((uint32x2_t*)a)[20]);
    t[1] = vmlal_high_u32(t[1], a[10], a[10]);
    t[0] = vmlal_u32(t[0], tmp[14], mod);
    t[1] = vmlal_u32(t[1], tmp[15], mod);
    out[12] = vmovn_u64( vandq_u64(t[0], mask));
    out[13] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], PER_LIMB);
    t[1] = vshrq_n_u64(t[1], PER_LIMB);
    //---x^21----
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[16], ((int32x2_t*)a)[26]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[8], ((int32x4_t*)a)[13]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[18], ((int32x2_t*)a)[24]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[9], ((int32x4_t*)a)[12]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[20], ((int32x2_t*)a)[22]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[10], ((int32x4_t*)a)[11]);
    t[0] = vmlal_u32(t[0], tmp[16], mod);
    t[1] = vmlal_u32(t[1], tmp[17], mod);
    out[14] = vmovn_u64( vandq_u64(t[0], mask));
    out[15] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], PER_LIMB);
    t[1] = vshrq_n_u64(t[1], PER_LIMB);
    //---x^22----
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[18], ((int32x2_t*)a)[26]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[9], ((int32x4_t*)a)[13]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[20], ((int32x2_t*)a)[24]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[10], ((int32x4_t*)a)[12]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[22], ((uint32x2_t*)a)[22]);
    t[1] = vmlal_high_u32(t[1], a[11], a[11]);
    t[0] = vmlal_u32(t[0], tmp[18], mod);
    t[1] = vmlal_u32(t[1], tmp[19], mod);
    out[16] = vmovn_u64( vandq_u64(t[0], mask));
    out[17] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], PER_LIMB);
    t[1] = vshrq_n_u64(t[1], PER_LIMB);
    //---x^23----
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[20], ((int32x2_t*)a)[26]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[10], ((int32x4_t*)a)[13]);
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[22], ((int32x2_t*)a)[24]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[11], ((int32x4_t*)a)[12]);
    t[0] = vmlal_u32(t[0], tmp[20], mod);
    t[1] = vmlal_u32(t[1], tmp[21], mod);
    out[18] = vmovn_u64( vandq_u64(t[0], mask));
    out[19] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], PER_LIMB);
    t[1] = vshrq_n_u64(t[1], PER_LIMB);
    //---x^24----
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[22], ((int32x2_t*)a)[26]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[11], ((int32x4_t*)a)[13]);
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[24], ((uint32x2_t*)a)[24]);
    t[1] = vmlal_high_u32(t[1], a[12], a[12]);
    t[0] = vmlal_u32(t[0], tmp[22], mod);
    t[1] = vmlal_u32(t[1], tmp[23], mod);
    out[20] = vmovn_u64( vandq_u64(t[0], mask));
    out[21] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], PER_LIMB);
    t[1] = vshrq_n_u64(t[1], PER_LIMB);
    //---x^25----
    t[0] = (uint64x2_t)vqdmlal_s32(((int64x2_t*)t)[0], ((int32x2_t*)a)[24], ((int32x2_t*)a)[26]);
    t[1] = (uint64x2_t)vqdmlal_high_s32(((int64x2_t*)t)[1], ((int32x4_t*)a)[12], ((int32x4_t*)a)[13]);
    t[0] = vmlal_u32(t[0], tmp[24], mod);
    t[1] = vmlal_u32(t[1], tmp[25], mod);
    out[22] = vmovn_u64( vandq_u64(t[0], mask));
    out[23] = vmovn_u64( vandq_u64(t[1], mask));
    t[0] = vshrq_n_u64(t[0], PER_LIMB);
    t[1] = vshrq_n_u64(t[1], PER_LIMB);
    //---x^26----
    t[0] = vmlal_u32(t[0], ((uint32x2_t*)a)[26], ((uint32x2_t*)a)[26]);
    t[1] = vmlal_high_u32(t[1], a[13], a[13]);
    t[0] = vmlal_u32(t[0], tmp[26], mod);
    t[1] = vmlal_u32(t[1], tmp[27], mod);
    out[24] = vmovn_u64( vandq_u64(t[0], mask));
    out[25] = vmovn_u64( vandq_u64(t[1], mask));


    //---x^27----
    out[26] = vshrn_n_u64(t[0], 28);
    out[27] = vshrn_n_u64(t[1], 28);
}

void fp_sqrt_batched(uint32x4_t *out, const uint32x4_t *in) {
    uint32x4_t x[FP_LIMBS];
    memmove(x, in, sizeof(uint32x4_t)*FP_LIMBS);
    fp_sqr_batched((uint32x2_t*)out, x);
    for (int i=0; i<373; i++){
      fp_sqr_batched((uint32x2_t*)out, out);
    }

    fp_sqr_batched((uint32x2_t*)x, out);
    for (int i=0; i<5; i++){
        fp_sqr_batched((uint32x2_t*)x, x);
    }

    //out = 2^374, x = 2^380
    fp_mul_batched((uint32x2_t*)out, out, x);
}

void fp_exp3div4_vec(uint32x4_t *out, const uint32x4_t *in) {
    uint32x4_t x[FP_LIMBS], t0[FP_LIMBS], t1[FP_LIMBS], t2[FP_LIMBS], t3[FP_LIMBS], t4[FP_LIMBS], t5[FP_LIMBS];
    memmove(x, in, sizeof(uint32x4_t)*FP_LIMBS);

    fp_sqr_batched((uint32x2_t *)out, x);      //out = in^2
    fp_sqr_batched((uint32x2_t *)t0, out);

    fp_mul_batched((uint32x2_t *)t1, x, t0);
    fp_mul_batched((uint32x2_t *)out, out, t1);
    fp_sqr_batched((uint32x2_t *)t0, out);
    fp_sqr_batched((uint32x2_t *)t3, t0);
    fp_sqr_batched((uint32x2_t *)t4, t3);
    fp_sqr_batched((uint32x2_t *)t2, t4);
    memmove(t5, t2, sizeof(uint32x4_t)*FP_LIMBS);
    for (int i=0; i<3; i++){
        fp_sqr_batched((uint32x2_t *)t5, t5);
    }
    fp_mul_batched((uint32x2_t *)t2, t2, t5);
    memmove(t5, t2, sizeof(uint32x4_t)*FP_LIMBS);
    for (int i=0; i<6; i++){
        fp_sqr_batched((uint32x2_t *)t5, t5);
    }
    fp_mul_batched((uint32x2_t *)t2, t2, t5);
    memmove(t5, t2, sizeof(uint32x4_t)*FP_LIMBS);
    for (int i=0; i<2; i++){
        fp_sqr_batched((uint32x2_t *)t5, t5);
    }
    fp_mul_batched((uint32x2_t *)t5, t4, t5);
    for (int i=0; i<13; i++){
        fp_sqr_batched((uint32x2_t *)t5, t5);
    }
    fp_mul_batched((uint32x2_t *)t2, t2, t5);
    memmove(t5, t2, sizeof(uint32x4_t)*FP_LIMBS);
    for (int i=0; i<2; i++){
        fp_sqr_batched((uint32x2_t *)t5, t5);
    }
    fp_mul_batched((uint32x2_t *)t4, t4, t5);
    for (int i=0; i<28; i++){
        fp_sqr_batched((uint32x2_t *)t4, t4);
    }
    fp_mul_batched((uint32x2_t *)t2, t2, t4);
    fp_sqr_batched((uint32x2_t *)t4, t2);
    fp_mul_batched((uint32x2_t *)t3, t3, t4);
    for (int i=0; i<59; i++){
        fp_sqr_batched((uint32x2_t *)t3, t3);
    }
    fp_mul_batched((uint32x2_t *)t2, t2, t3);
    fp_mul_batched((uint32x2_t *)t1, t1, t2);
    fp_mul_batched((uint32x2_t *)out, out, t1);
    fp_mul_batched((uint32x2_t *)t0, t0, out);
    fp_mul_batched((uint32x2_t *)t1, t0, t1);
    fp_sqr_batched((uint32x2_t *)t2, t1);
    fp_mul_batched((uint32x2_t *)t2, t1, t2);
    fp_sqr_batched((uint32x2_t *)t2, t2);
    fp_mul_batched((uint32x2_t *)t2, t1, t2);
    fp_mul_batched((uint32x2_t *)t0, t0, t2);
    fp_mul_batched((uint32x2_t *)out, out, t0);
    fp_sqr_batched((uint32x2_t *)t2, out);
    fp_mul_batched((uint32x2_t *)t2, out, t2);
    fp_mul_batched((uint32x2_t *)t0, t0, t2);
    fp_mul_batched((uint32x2_t *)t1, t1, t0);
    memmove(t2, t1, sizeof(uint32x4_t)*FP_LIMBS);
    for (int i=0; i<128; i++){
        fp_sqr_batched((uint32x2_t *)t2, t2);
    }
    fp_mul_batched((uint32x2_t *)t1, t1, t2);
    fp_mul_batched((uint32x2_t *)t0, t0, t1);
    for (int i=0; i<125; i++){
        fp_sqr_batched((uint32x2_t *)t0, t0);
    }
    fp_mul_batched((uint32x2_t *)out, out, t0);
}

void fp_neg_vec(uint32x4_t *out, const uint32x4_t *in){
    uint32x4_t x[FP_LIMBS], q[FP_LIMBS];
    memmove(x, in, sizeof(uint32x4_t)*FP_LIMBS);
    for(int i = 0;i<(FP_LIMBS-1);i++) q[i] = vdupq_n_u32(Q2_VALUE);
    q[FP_LIMBS-1] = vdupq_n_u32(Q2_VALUE_SIGNIFICANT);
    
    for (int i=0; i<FP_LIMBS; i++){
        out[i] = vsubq_u32(q[i], x[i]);
    }
    fp_bactched_reduction(out);
}

void modmul32(const uint32_t *a, const uint32_t *b, uint32_t *c) {
    uint64_t t = 0;
    uint32_t p8 = MASK_VALUE;
    uint32_t q = ((uint32_t)1 << 28u); // q is unsaturated radix
    uint32_t mask = (uint32_t)(q - (uint32_t)1);
    
    //---x^0----
    t += (uint64_t)a[0] * b[0];
    uint32_t v0 = ((uint32_t)t & mask);
    t >>= PER_LIMB;
    //---x^1----
    t += (uint64_t)a[0] * b[1];
    t += (uint64_t)a[1] * b[0];
    uint32_t v1 = ((uint32_t)t & mask);
    t >>= PER_LIMB;
    //---x^2----
    t += (uint64_t)a[0] * b[2];
    t += (uint64_t)a[1] * b[1];
    t += (uint64_t)a[2] * b[0];
    uint32_t v2 = ((uint32_t)t & mask);
    t >>= PER_LIMB;
    //---x^3----
    t += (uint64_t)a[0] * b[3];
    t += (uint64_t)a[1] * b[2];
    t += (uint64_t)a[2] * b[1];
    t += (uint64_t)a[3] * b[0];
    uint32_t v3 = ((uint32_t)t & mask);
    t >>= PER_LIMB;
    //---x^4----
    t += (uint64_t)a[0] * b[4];
    t += (uint64_t)a[1] * b[3];
    t += (uint64_t)a[2] * b[2];
    t += (uint64_t)a[3] * b[1];
    t += (uint64_t)a[4] * b[0];
    uint32_t v4 = ((uint32_t)t & mask);
    t >>= PER_LIMB;
    //---x^5----
    t += (uint64_t)a[0] * b[5];
    t += (uint64_t)a[1] * b[4];
    t += (uint64_t)a[2] * b[3];
    t += (uint64_t)a[3] * b[2];
    t += (uint64_t)a[4] * b[1];
    t += (uint64_t)a[5] * b[0];
    uint32_t v5 = ((uint32_t)t & mask);
    t >>= PER_LIMB;
    //---x^6----
    t += (uint64_t)a[0] * b[6];
    t += (uint64_t)a[1] * b[5];
    t += (uint64_t)a[2] * b[4];
    t += (uint64_t)a[3] * b[3];
    t += (uint64_t)a[4] * b[2];
    t += (uint64_t)a[5] * b[1];
    t += (uint64_t)a[6] * b[0];
    uint32_t v6 = ((uint32_t)t & mask);
    t >>= PER_LIMB;
    //---x^7----
    t += (uint64_t)a[0] * b[7];
    t += (uint64_t)a[1] * b[6];
    t += (uint64_t)a[2] * b[5];
    t += (uint64_t)a[3] * b[4];
    t += (uint64_t)a[4] * b[3];
    t += (uint64_t)a[5] * b[2];
    t += (uint64_t)a[6] * b[1];
    t += (uint64_t)a[7] * b[0];
    uint32_t v7 = ((uint32_t)t & mask);
    t >>= PER_LIMB;
    //---x^8----
    t += (uint64_t)a[0] * b[8];
    t += (uint64_t)a[1] * b[7];
    t += (uint64_t)a[2] * b[6];
    t += (uint64_t)a[3] * b[5];
    t += (uint64_t)a[4] * b[4];
    t += (uint64_t)a[5] * b[3];
    t += (uint64_t)a[6] * b[2];
    t += (uint64_t)a[7] * b[1];
    t += (uint64_t)a[8] * b[0];
    uint32_t v8 = ((uint32_t)t & mask);
    t >>= PER_LIMB;
    //---x^9----
    t += (uint64_t)a[0] * b[9];
    t += (uint64_t)a[1] * b[8];
    t += (uint64_t)a[2] * b[7];
    t += (uint64_t)a[3] * b[6];
    t += (uint64_t)a[4] * b[5];
    t += (uint64_t)a[5] * b[4];
    t += (uint64_t)a[6] * b[3];
    t += (uint64_t)a[7] * b[2];
    t += (uint64_t)a[8] * b[1];
    t += (uint64_t)a[9] * b[0];
    uint32_t v9 = ((uint32_t)t & mask);
    t >>= PER_LIMB;
    //---x^10----
    t += (uint64_t)a[0] * b[10];
    t += (uint64_t)a[1] * b[9];
    t += (uint64_t)a[2] * b[8];
    t += (uint64_t)a[3] * b[7];
    t += (uint64_t)a[4] * b[6];
    t += (uint64_t)a[5] * b[5];
    t += (uint64_t)a[6] * b[4];
    t += (uint64_t)a[7] * b[3];
    t += (uint64_t)a[8] * b[2];
    t += (uint64_t)a[9] * b[1];
    t += (uint64_t)a[10] * b[0];
    uint32_t v10 = ((uint32_t)t & mask);
    t >>= PER_LIMB;
    //---x^11----
    t += (uint64_t)a[0] * b[11];
    t += (uint64_t)a[1] * b[10];
    t += (uint64_t)a[2] * b[9];
    t += (uint64_t)a[3] * b[8];
    t += (uint64_t)a[4] * b[7];
    t += (uint64_t)a[5] * b[6];
    t += (uint64_t)a[6] * b[5];
    t += (uint64_t)a[7] * b[4];
    t += (uint64_t)a[8] * b[3];
    t += (uint64_t)a[9] * b[2];
    t += (uint64_t)a[10] * b[1];
    t += (uint64_t)a[11] * b[0];
    uint32_t v11 = ((uint32_t)t & mask);
    t >>= PER_LIMB;
    //---x^12----
    t += (uint64_t)a[0] * b[12];
    t += (uint64_t)a[1] * b[11];
    t += (uint64_t)a[2] * b[10];
    t += (uint64_t)a[3] * b[9];
    t += (uint64_t)a[4] * b[8];
    t += (uint64_t)a[5] * b[7];
    t += (uint64_t)a[6] * b[6];
    t += (uint64_t)a[7] * b[5];
    t += (uint64_t)a[8] * b[4];
    t += (uint64_t)a[9] * b[3];
    t += (uint64_t)a[10] * b[2];
    t += (uint64_t)a[11] * b[1];
    t += (uint64_t)a[12] * b[0];
    uint32_t v12 = ((uint32_t)t & mask);
    t >>= PER_LIMB;
    //---x^13----
    t += (uint64_t)a[0] * b[13];
    t += (uint64_t)a[1] * b[12];
    t += (uint64_t)a[2] * b[11];
    t += (uint64_t)a[3] * b[10];
    t += (uint64_t)a[4] * b[9];
    t += (uint64_t)a[5] * b[8];
    t += (uint64_t)a[6] * b[7];
    t += (uint64_t)a[7] * b[6];
    t += (uint64_t)a[8] * b[5];
    t += (uint64_t)a[9] * b[4];
    t += (uint64_t)a[10] * b[3];
    t += (uint64_t)a[11] * b[2];
    t += (uint64_t)a[12] * b[1];
    t += (uint64_t)a[13] * b[0];
    t += (uint64_t)v0 * (uint64_t)p8;
    uint32_t v13 = ((uint32_t)t & mask);
    t >>= PER_LIMB;
    //---x^14----
    t += (uint64_t)a[1] * b[13];
    t += (uint64_t)a[2] * b[12];
    t += (uint64_t)a[3] * b[11];
    t += (uint64_t)a[4] * b[10];
    t += (uint64_t)a[5] * b[9];
    t += (uint64_t)a[6] * b[8];
    t += (uint64_t)a[7] * b[7];
    t += (uint64_t)a[8] * b[6];
    t += (uint64_t)a[9] * b[5];
    t += (uint64_t)a[10] * b[4];
    t += (uint64_t)a[11] * b[3];
    t += (uint64_t)a[12] * b[2];
    t += (uint64_t)a[13] * b[1];
    t += (uint64_t)v1 * (uint64_t)p8;
    c[0] = ((uint32_t)t & mask);
    t >>= PER_LIMB;
    //---x^15----
    t += (uint64_t)a[2] * b[13];
    t += (uint64_t)a[3] * b[12];
    t += (uint64_t)a[4] * b[11];
    t += (uint64_t)a[5] * b[10];
    t += (uint64_t)a[6] * b[9];
    t += (uint64_t)a[7] * b[8];
    t += (uint64_t)a[8] * b[7];
    t += (uint64_t)a[9] * b[6];
    t += (uint64_t)a[10] * b[5];
    t += (uint64_t)a[11] * b[4];
    t += (uint64_t)a[12] * b[3];
    t += (uint64_t)a[13] * b[2];
    t += (uint64_t)v2 * (uint64_t)p8;
    c[1] = ((uint32_t)t & mask);
    t >>= PER_LIMB;
    //---x^16----
    t += (uint64_t)a[3] * b[13];
    t += (uint64_t)a[4] * b[12];
    t += (uint64_t)a[5] * b[11];
    t += (uint64_t)a[6] * b[10];
    t += (uint64_t)a[7] * b[9];
    t += (uint64_t)a[8] * b[8];
    t += (uint64_t)a[9] * b[7];
    t += (uint64_t)a[10] * b[6];
    t += (uint64_t)a[11] * b[5];
    t += (uint64_t)a[12] * b[4];
    t += (uint64_t)a[13] * b[3];
    t += (uint64_t)v3 * (uint64_t)p8;
    c[2] = ((uint32_t)t & mask);
    t >>= PER_LIMB;
    //---x^17----
    t += (uint64_t)a[4] * b[13];
    t += (uint64_t)a[5] * b[12];
    t += (uint64_t)a[6] * b[11];
    t += (uint64_t)a[7] * b[10];
    t += (uint64_t)a[8] * b[9];
    t += (uint64_t)a[9] * b[8];
    t += (uint64_t)a[10] * b[7];
    t += (uint64_t)a[11] * b[6];
    t += (uint64_t)a[12] * b[5];
    t += (uint64_t)a[13] * b[4];
    t += (uint64_t)v4 * (uint64_t)p8;
    c[3] = ((uint32_t)t & mask);
    t >>= PER_LIMB;
    //---x^18----
    t += (uint64_t)a[5] * b[13];
    t += (uint64_t)a[6] * b[12];
    t += (uint64_t)a[7] * b[11];
    t += (uint64_t)a[8] * b[10];
    t += (uint64_t)a[9] * b[9];
    t += (uint64_t)a[10] * b[8];
    t += (uint64_t)a[11] * b[7];
    t += (uint64_t)a[12] * b[6];
    t += (uint64_t)a[13] * b[5];
    t += (uint64_t)v5 * (uint64_t)p8;
    c[4] = ((uint32_t)t & mask);
    t >>= PER_LIMB;
    //---x^19----
    t += (uint64_t)a[6] * b[13];
    t += (uint64_t)a[7] * b[12];
    t += (uint64_t)a[8] * b[11];
    t += (uint64_t)a[9] * b[10];
    t += (uint64_t)a[10] * b[9];
    t += (uint64_t)a[11] * b[8];
    t += (uint64_t)a[12] * b[7];
    t += (uint64_t)a[13] * b[6];
    t += (uint64_t)v6 * (uint64_t)p8;
    c[5] = ((uint32_t)t & mask);
    t >>= PER_LIMB;
    //---x^20----
    t += (uint64_t)a[7] * b[13];
    t += (uint64_t)a[8] * b[12];
    t += (uint64_t)a[9] * b[11];
    t += (uint64_t)a[10] * b[10];
    t += (uint64_t)a[11] * b[9];
    t += (uint64_t)a[12] * b[8];
    t += (uint64_t)a[13] * b[7];
    t += (uint64_t)v7 * (uint64_t)p8;
    c[6] = ((uint32_t)t & mask);
    t >>= PER_LIMB;
    //---x^21----
    t += (uint64_t)a[8] * b[13];
    t += (uint64_t)a[9] * b[12];
    t += (uint64_t)a[10] * b[11];
    t += (uint64_t)a[11] * b[10];
    t += (uint64_t)a[12] * b[9];
    t += (uint64_t)a[13] * b[8];
    t += (uint64_t)v8 * (uint64_t)p8;
    c[7] = ((uint32_t)t & mask);
    t >>= PER_LIMB;
    //---x^22----
    t += (uint64_t)a[9] * b[13];
    t += (uint64_t)a[10] * b[12];
    t += (uint64_t)a[11] * b[11];
    t += (uint64_t)a[12] * b[10];
    t += (uint64_t)a[13] * b[9];
    t += (uint64_t)v9 * (uint64_t)p8;
    c[8] = ((uint32_t)t & mask);
    t >>= PER_LIMB;
    //---x^23----
    t += (uint64_t)a[10] * b[13];
    t += (uint64_t)a[11] * b[12];
    t += (uint64_t)a[12] * b[11];
    t += (uint64_t)a[13] * b[10];
    t += (uint64_t)v10 * (uint64_t)p8;
    c[9] = ((uint32_t)t & mask);
    t >>= PER_LIMB;
    //---x^24----
    t += (uint64_t)a[11] * b[13];
    t += (uint64_t)a[12] * b[12];
    t += (uint64_t)a[13] * b[11];
    t += (uint64_t)v11 * (uint64_t)p8;
    c[10] = ((uint32_t)t & mask);
    t >>= PER_LIMB;
    //---x^25----
    t += (uint64_t)a[12] * b[13];
    t += (uint64_t)a[13] * b[12];
    t += (uint64_t)v12 * (uint64_t)p8;
    c[11] = ((uint32_t)t & mask);
    t >>= PER_LIMB;
    //---x^26----
    t += (uint64_t)a[13] * b[13];
    t += (uint64_t)v13 * (uint64_t)p8;
    c[12] = ((uint32_t)t & mask);
    t >>= PER_LIMB;
    c[13] = (uint32_t)t;
}

uint32_t prop32 (uint32_t *n) {
    int i;
    uint32_t mask = ((uint32_t)1 << 28u) - (uint32_t)1;
    int32_t carry = (int32_t)n[0];
    carry >>= 28u;
    n[0] &= mask;
    for (i = 1; i < (FP_LIMBS-1); i++) {
        carry += (int32_t)n[i];
        n[i] = (uint32_t)carry & mask;
        carry >>= 28u;
    }
    n[FP_LIMBS-1] += (uint32_t)carry;
    return -((n[FP_LIMBS-1] >> 1) >> 30u);
}

int flatten32(uint32_t *n) {
    uint32_t carry = prop32(n);
    n[0] -= (uint32_t)1u & carry;
    n[FP_LIMBS-1] += ((uint32_t)MASK_VALUE) & carry;
    (void)prop32(n);
    return (int)(carry & 1);
}

int modfsb32(uint32_t *n) {
  n[0] += (uint32_t)1u;
  n[FP_LIMBS-1] -= (uint32_t)MASK_VALUE;
  return flatten32(n);
}

void redc32(uint32_t *n, uint32_t *m) {
  int i;
  uint32_t c[FP_LIMBS];
  c[0] = 1;
  for (i = 1; i < FP_LIMBS; i++) {
    c[i] = 0;
  }
  modmul32(n, c, m);
  (void)modfsb32(m);
}

uint32_t fp_is_zero_32(uint32_t* p){
  uint32_t c[FP_LIMBS], d = 0;
  redc32(p, c);
  for (int i = 0; i < FP_LIMBS; i++) {
    d |= c[i];
  }
  return -(uint32_t)((uint32_t)1 & ((d - (uint32_t)1) >> 28u));
}

uint32x4_t theta_point_is_zero(const uint32x4_t* a){
    uint32x4_t mask = (uint32x4_t)vdupq_n_s32(-1);
    uint32x4_t zero = vdupq_n_u32(0);
    uint32x4_t qlow = vdupq_n_u32((1<<PER_LIMB)-1);
    uint32x4_t qhigh = vdupq_n_u32(MASK_VALUE-1);
    uint32x4_t tmp;
    
    for(int i = 0;i<(FP_LIMBS-1);i++){
        tmp = vorrq_u32(vceqq_u32(a[i], qlow), vceqq_u32(a[i], zero));
        mask = vandq_u32(mask, tmp);
    }
    tmp = vorrq_u32(vceqq_u32(a[FP_LIMBS-1], qhigh), vceqq_u32(a[FP_LIMBS-1], zero));
    return vandq_u32(mask, tmp);
}