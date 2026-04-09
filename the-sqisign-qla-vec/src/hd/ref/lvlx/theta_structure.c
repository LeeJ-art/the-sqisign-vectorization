#include "theta_structure.h"
#include <assert.h>
#include <arm_neon.h>

void
theta_precomputation(theta_structure_t *A)
{

    if (A->precomputation) {
        return;
    }

    theta_point_t A_dual;
    to_squared_theta(&A_dual, &A->null_point);

    fp2_t t1, t2;
    fp2_mul(&t1, &A_dual.x, &A_dual.y);
    fp2_mul(&t2, &A_dual.z, &A_dual.t);
    fp2_mul(&A->XYZ0, &t1, &A_dual.z);
    fp2_mul(&A->XYT0, &t1, &A_dual.t);
    fp2_mul(&A->YZT0, &t2, &A_dual.y);
    fp2_mul(&A->XZT0, &t2, &A_dual.x);

    fp2_mul(&t1, &A->null_point.x, &A->null_point.y);
    fp2_mul(&t2, &A->null_point.z, &A->null_point.t);
    fp2_mul(&A->xyz0, &t1, &A->null_point.z);
    fp2_mul(&A->xyt0, &t1, &A->null_point.t);
    fp2_mul(&A->yzt0, &t2, &A->null_point.y);
    fp2_mul(&A->xzt0, &t2, &A->null_point.x);

    A->precomputation = true;
}

void structure2point(theta_point_t* out, theta_structure_t *A){
    out[0].x = A->XYZ0;
    out[0].y = A->YZT0;
    out[0].z = A->XZT0;
    out[0].t = A->XYT0;
    out[1].x = A->xyz0;
    out[1].y = A->yzt0;
    out[1].z = A->xzt0;
    out[1].t = A->xyt0;
}

void structure2point_reindex(theta_point_t* out, theta_structure_t *A){
    out[0].x = A->YZT0;
    out[0].y = A->XZT0;
    out[0].z = A->XYT0;
    out[0].t = A->XYZ0;
    out[1].x = A->yzt0;
    out[1].y = A->xzt0;
    out[1].z = A->xyt0;
    out[1].t = A->xyz0;
}

void point2structure(theta_structure_t* A, theta_point_t *out){
    A->YZT0 = out[0].x;
    A->XZT0 = out[0].y;
    A->XYT0 = out[0].z;
    A->XYZ0 = out[0].t;

    A->yzt0 = out[1].x;
    A->xzt0 = out[1].y;
    A->xyt0 = out[1].z;
    A->xyz0 = out[1].t;
}

void theta_precomputation_vec(uint32x4_t* a0, uint32x4_t* a1, uint32x4_t* A_null, bool* A_precomputeFlag)
{
    if (A_precomputeFlag[0]) {
        return;
    }
    
    // theta_point_t A_dual;
    uint32x4_t A_dual[FP2_LIMBS];
    // to_squared_theta(&A_dual, &A->null_point);
    to_squared_theta_batched(A_dual, A_null);

    // fp2_t t1, t2;
    uint32x4_t t1[FP2_LIMBS], t2[FP2_LIMBS], t3[FP2_LIMBS];
    // reindex
    for(int i = 0;i<FP2_LIMBS;i++){
        t1[i][0] = A_dual[i][0];
        t1[i][1] = A_dual[i][2];
        t1[i][2] = A_null[i][0];
        t1[i][3] = A_null[i][2];
        t2[i][0] = A_dual[i][1];
        t2[i][1] = A_dual[i][3];
        t2[i][2] = A_null[i][1];
        t2[i][3] = A_null[i][3];
    }
    // t3 = XY ZT xy zt
    fp2_mul_batched(t3, t1, t2);
    for(int i = 0;i<FP2_LIMBS;i++){
        t1[i][0] = t3[i][1];
        t1[i][1] = t3[i][1];
        t1[i][2] = t3[i][0];
        t1[i][3] = t3[i][0];
        t2[i][0] = A_dual[i][1];
        t2[i][1] = A_dual[i][0];
        t2[i][2] = A_dual[i][3];
        t2[i][3] = A_dual[i][2];
    }
    // a0 = {YZT, XZT, XYT, XYZ}
    fp2_mul_batched(a0, t1, t2);
    
    for(int i = 0;i<FP2_LIMBS;i++){
        t1[i][0] = t3[i][3];
        t1[i][1] = t3[i][3];
        t1[i][2] = t3[i][2];
        t1[i][3] = t3[i][2];
        t2[i][0] = A_null[i][1];
        t2[i][1] = A_null[i][0];
        t2[i][2] = A_null[i][3];
        t2[i][3] = A_null[i][2];
    }
    // a1 = {yzt, xzt, xyt, xyz}
    fp2_mul_batched(a1, t1, t2);

    memcpy(t1, a0, sizeof(uint32x4_t)*FP2_LIMBS);
    memcpy(t2, a1, sizeof(uint32x4_t)*FP2_LIMBS);
    for (int i=0; i<FP2_LIMBS; i++){
        a0[i][0] = t1[i][3];
        a0[i][1] = t1[i][0];
        a0[i][2] = t1[i][1];
        a0[i][3] = t1[i][2];

        a1[i][0] = t2[i][3];
        a1[i][1] = t2[i][0];
        a1[i][2] = t2[i][1];
        a1[i][3] = t2[i][2];
    }

    A_precomputeFlag[0] = true;
}



void
double_point(theta_point_t *out, theta_structure_t *A, const theta_point_t *in)
{
    to_squared_theta(out, in);
    fp2_sqr(&out->x, &out->x);
    fp2_sqr(&out->y, &out->y);
    fp2_sqr(&out->z, &out->z);
    fp2_sqr(&out->t, &out->t);

    if (!A->precomputation) {
        theta_precomputation(A);
    }
    fp2_mul(&out->x, &out->x, &A->YZT0);
    fp2_mul(&out->y, &out->y, &A->XZT0);
    fp2_mul(&out->z, &out->z, &A->XYT0);
    fp2_mul(&out->t, &out->t, &A->XYZ0);

    hadamard(out, out);

    fp2_mul(&out->x, &out->x, &A->yzt0);
    fp2_mul(&out->y, &out->y, &A->xzt0);
    fp2_mul(&out->z, &out->z, &A->xyt0);
    fp2_mul(&out->t, &out->t, &A->xyz0);
}

void double_point_vec(uint32x4_t *out, 
    uint32x4_t *A_null_point, uint32x4_t *A_capital, uint32x4_t *A_lowercase, bool *A_precomputeFlag, 
    const uint32x4_t *in)
{
    //memcpy
    uint32x4_t p1[FP2_LIMBS];
    memcpy(p1, in, sizeof(uint32x4_t)*FP2_LIMBS);
    
    // to_squared_theta(out, in);
    to_squared_theta_batched(out, p1);

    fp2_sqr_batched(out, out);

    // Here keep A's points in a0, a1 with mul-friendly sequence{YZT, XZT, XYT, XYZ}
    if (!A_precomputeFlag[0]) {
        // theta_precomputation(A);
        theta_precomputation_vec(A_capital, A_lowercase, A_null_point, A_precomputeFlag);
    }

    uint32x4_t a0[FP2_LIMBS], a1[FP2_LIMBS];
    //re-index
    for (int i=0; i<FP2_LIMBS; i++){
        a0[i][0] = A_capital[i][1];
        a0[i][1] = A_capital[i][2];
        a0[i][2] = A_capital[i][3];
        a0[i][3] = A_capital[i][0];

        a1[i][0] = A_lowercase[i][1];
        a1[i][1] = A_lowercase[i][2];
        a1[i][2] = A_lowercase[i][3];
        a1[i][3] = A_lowercase[i][0];
    }

    // fp2_mul(&out->x, &out->x, &A->YZT0);
    // fp2_mul(&out->y, &out->y, &A->XZT0);
    // fp2_mul(&out->z, &out->z, &A->XYT0);
    // fp2_mul(&out->t, &out->t, &A->XYZ0);
    fp2_mul_batched(out, out, a0);

    //hadamard(out, out);
    hadamard_vec(out, out);

    // fp2_mul(&out->x, &out->x, &A->yzt0);
    // fp2_mul(&out->y, &out->y, &A->xzt0);
    // fp2_mul(&out->z, &out->z, &A->xyt0);
    // fp2_mul(&out->t, &out->t, &A->xyz0);
    fp2_bactched_reduction(out);
    fp2_bactched_reduction(a1);
    fp2_mul_batched(out, out, a1);
}



void
double_iter(theta_point_t *out, theta_structure_t *A, const theta_point_t *in, int exp)
{
    if (exp == 0) {
        *out = *in;
    } else {
        double_point(out, A, in);
        for (int i = 1; i < exp; i++) {
            double_point(out, A, out);
        }
    }
}

void double_iter_vec(uint32x4_t *out, 
     uint32x4_t *theta_structure_null_point, uint32x4_t *theta_structure_capital, uint32x4_t *theta_structure_lowercase, bool *theta_structure_precomputeFlag, 
    const uint32x4_t *in, const int exp)
{
    if (exp == 0) {
        memcpy(out, in, sizeof(uint32x4_t)*FP2_LIMBS);        
    } else {
        double_point_vec(out, theta_structure_null_point, theta_structure_capital, theta_structure_lowercase, theta_structure_precomputeFlag, in);
        for (int i = 1; i < exp; i++) {
            double_point_vec(out, theta_structure_null_point, theta_structure_capital, theta_structure_lowercase, theta_structure_precomputeFlag, out);
        }
    }
}

uint32_t
is_product_theta_point(const theta_point_t *P)
{
    fp2_t t1, t2;
    fp2_mul(&t1, &P->x, &P->t);
    fp2_mul(&t2, &P->y, &P->z);
    return fp2_is_equal(&t1, &t2);
}

uint32_t 
is_product_theta_point_vec(const uint32x4_t *P)
{
    // Returns 0xFF...FF (true) if a=b, 0 (false) otherwise
    uint32x4_t a[FP2_LIMBS], b[FP2_LIMBS];
    uint32_t flag = 1;
    memmove(a, P, sizeof(uint32x4_t)*FP2_LIMBS);
    for (int i=0; i<FP2_LIMBS; i++){
        b[i][0] = a[i][3];
        b[i][1] = a[i][2];
    }
    // fp2_mul(&t1, &P->x, &P->t);
    // fp2_mul(&t2, &P->y, &P->z);
    fp2_mul_batched(a, a, b);
    for (int i=0; i<FP2_LIMBS; i++){
        flag &= (((a[i][0] ^ a[i][1]) - 1) >> PER_LIMB) & 1;
    }
    return -(flag);
}