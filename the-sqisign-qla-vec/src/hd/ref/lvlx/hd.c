#include <hd.h>
#include <mp.h>
#include <assert.h>
#include <stdio.h>
#include <bench.h>
#include <arm_neon.h>

static __inline__ uint64_t rdtsc(void)
{
    return (uint64_t)cpucycles();
}


void
double_couple_point(theta_couple_point_t *out, const theta_couple_point_t *in, const theta_couple_curve_t *E1E2)
{
    ec_dbl(&out->P1, &in->P1, &E1E2->E1);
    ec_dbl(&out->P2, &in->P2, &E1E2->E2);
}

void
double_couple_point_iter(theta_couple_point_t *out,
                         unsigned n,
                         const theta_couple_point_t *in,
                         const theta_couple_curve_t *E1E2)
{
    if (n == 0) {
        memmove(out, in, sizeof(theta_couple_point_t));
    } else {
        double_couple_point(out, in, E1E2);
        for (unsigned i = 0; i < n - 1; i++) {
            double_couple_point(out, out, E1E2);
        }
    }
}

void
add_couple_jac_points(theta_couple_jac_point_t *out,
                      const theta_couple_jac_point_t *T1,
                      const theta_couple_jac_point_t *T2,
                      const theta_couple_curve_t *E1E2)
{
    ADD(&out->P1, &T1->P1, &T2->P1, &E1E2->E1);
    ADD(&out->P2, &T1->P2, &T2->P2, &E1E2->E2);
}

void
double_couple_jac_point(theta_couple_jac_point_t *out,
                        const theta_couple_jac_point_t *in,
                        const theta_couple_curve_t *E1E2)
{
    DBL(&out->P1, &in->P1, &E1E2->E1);
    DBL(&out->P2, &in->P2, &E1E2->E2);
}

void DBL_vec(uint32x4_t *out1, uint32x4_t *out2, const uint32x4_t *jac1, const uint32x4_t *jac2,
                        const uint32x4_t *E1, const uint32x4_t *E2)
{ // Cost of 6M + 6S.
  // Doubling on a Montgomery curve, representation in Jacobian coordinates (X:Y:Z) corresponding to
  // (X/Z^2,Y/Z^3) This version receives the coefficient value A

    /*P1.x with r=2**261*/
    uint32x4_t ax[FP2_LIMBS], az[FP2_LIMBS], a0[FP2_LIMBS], a1[FP2_LIMBS], a2[FP2_LIMBS], a3[FP2_LIMBS], a4[FP2_LIMBS], a5[FP2_LIMBS];
    for (int i=0; i<FP2_LIMBS; i++){      //ax = [x1, x2, z1, z2]
        ax[i][0] = jac1[i][0];
        ax[i][1] = jac2[i][0];
        ax[i][2] = jac1[i][2];
        ax[i][3] = jac2[i][2];
    }
    fp2_sqr_batched(a0, ax); //a0 = [x^2, x^2, z^2, z^2]


    for (int i=0; i<FP2_LIMBS; i++){      //az = [y1, y2, A1, A2]
        az[i][0] = jac1[i][1];
        az[i][1] = jac2[i][1];
        az[i][2] = E1[i][0];
        az[i][3] = E2[i][0];
    }
    for (int i=0; i<FP2_LIMBS; i++) a1[i] = vextq_u32(a0[i], az[i], 2);
    fp2_sqr_batched(a1, a1); //a1 = [z^4, z^4, y^2, y^2]
         
    for (int i=0; i<FP2_LIMBS; i++){
        for (int j=0; j<2; j++){
            a2[i][j] = a0[i][j];
            a2[i][j+2] = a1[i][j+2];
        }
    }
    fp2_add_batched(a1, a1, a2);         // a1 = [alpha, alpha, 2y^2, 2y^2]
    fp2_bactched_reduction(a1);

    for (int i=0; i<FP2_LIMBS; i++){
        for (int j=0; j<2; j++){
            a2[i][j+2] = ax[i][j];
            a0[i][j]   = a0[i][j+2];     // a0 = [z^2, z^2, 2y^2, 2y^2]
            a0[i][j+2] = a1[i][j+2]; 
        }
    }
    fp2_add_batched(a2, a2, a2);         // a2 = [2x^2, 2x^2, 2x, 2x]

    for (int i=0; i<FP2_LIMBS; i++) a3[i] = vextq_u32(az[i], ax[i], 2);
    fp2_mul_batched(a0, a0, a3);         // a0 = [Az^2, Az^2, 2y^2x, 2y^2x]


    for (int i=0; i<FP2_LIMBS; i++) a4[i] = vextq_u32(a0[i], a2[i], 2);
    for (int i=0; i<FP2_LIMBS; i++) a5[i] = vextq_u32(a0[i], a1[i], 2);
    fp2_add_batched(a3, a4, a5);         // a3 = [4y^2x, 4y^2x, 3x^2+z^4, 3x^2+z^4]

    for (int i=0; i<FP2_LIMBS; i++){
        for (int j=0; j<2; j++){
            a0[i][j+2] = a5[i][j+2] = a3[i][j];  //a0 = [Az^2, Az^2, 4y^2x, 4y^2x]
            a5[i][j] = a2[i][j+2];
        }
    }
    fp2_add_batched(a2, a0, a5);         // a2 = [(Az1^2 + 2x), (Az1^2 + 2x), 8y^2x, 8y^2x]

    for (int i=0; i<FP2_LIMBS; i++) a5[i] = vextq_u32(a1[i], a5[i], 2); //a5 = [2y^2, 2y^2, 2x, 2x]
    for (int i=0; i<FP2_LIMBS; i++){
        for (int j=0; j<2; j++){
            a4[i][j] = a2[i][j];  //a4 = [(Az1^2 + 2x), (Az1^2 + 2x), Az^2, Az^2]
            a4[i][j+2] = a0[i][j];
        }
    }
    fp2_bactched_reduction(a5);
    fp2_bactched_reduction(a4);
    fp2_mul_batched(a4, a4, a5);         // a4 = [2y^2(Az1^2 + 2x), 2y^2(Az1^2 + 2x), 2xAz^2, 2xAz^2]

    for (int i=0; i<FP2_LIMBS; i++){
        for (int j=0; j<2; j++){
            a5[i][j] = a4[i][j];
            a5[i][j+2] = a3[i][j+2];  //a5 = [2y^2(Az1^2 + 2x), 2y^2(Az1^2 + 2x), 3x^2+z^4, 3x^2+z^4]
        }
    }
    fp2_add_batched(a4, a4, a5);      // a4 = [4y^2(Az1^2 + 2x), 4y^2(Az1^2 + 2x), alpha, alpha]
    fp2_bactched_reduction(a4);

    for (int i=0; i<FP2_LIMBS; i++){
        for (int j=0; j<2; j++){
            a0[i][j] = a4[i][j+2];
            a0[i][j+2] = a1[i][j+2];
        }
    }
    fp2_sqr_batched(a0, a0);  //a0 = [alpha^2, alpha^2, 4y^4, 4y^4]

    uint32_t q[FP_LIMBS];
    for (int i=0; i<(FP_LIMBS-1); i++) q[i] = Q1_VALUE;
    q[FP_LIMBS-1] = Q1_VALUE_SIGNIFICANT;
    for (int i=0; i<FP2_LIMBS; i++){
        for (int j=0; j<2; j++){
            a1[i][j] = a4[i][j];
            a1[i][j+2] = a0[i][j+2];     //a1 = [4y^2(Az1^2 + 2x), 4y^2(Az1^2 + 2x), 4y^4, 4y^4]
            a0[i][j+2] = q[i%FP_LIMBS];
        }
    }
    fp2_sub_batched(a0, a0, a1);        // a0 = [alpha^2 - 4y^2(Az1^2 + 2x), alpha^2 - 4y^2(Az1^2 + 2x), -4y1^4, -4y1^4]
    fp2_bactched_reduction(a0);

    for (int i=0; i<FP2_LIMBS; i++){
        for (int j=0; j<2; j++){
            a3[i][j+2] = a0[i][j+2];     //a3 = [4y^2x, 4y^2x, -4y1^4, -4y1^4]
            a0[i][j+2] = a1[i][j+2];
        }
    }
    fp2_sub_batched(a3, a3, a0);        //a3 = [(4y^2x - Qx), (4y^2x - Qx), -8y1^4, -8y1^4]
    fp2_bactched_reduction(a3);

    for (int i=0; i<FP2_LIMBS; i++) a5[i] = vextq_u32(a4[i], az[i], 2); //a5 = [alpha, alpha, y, y]
    for (int i=0; i<FP2_LIMBS; i++){
        for (int j=0; j<2; j++){
            ax[i][j] = a3[i][j];     //ax = [(4y^2x - Qx), (4y^2x - Qx), z, z]     
        }
    }
    fp2_mul_batched(a1, a5, ax);         // a1 = [Qy, Qy, yz, yz]

    for (int i=0; i<FP2_LIMBS; i++){
        for (int j=0; j<2; j++){
            ax[i][j] = a3[i][j+2];     //ax = [-8y1^4, -8y1^4, yz, yz]
            ax[i][j+2] = a1[i][j+2];      
        }
    }
    fp2_add_batched(a1, a1, ax);        //a1 = [Qy, Qy, Qz, Qz]
    fp2_bactched_reduction(a1);


    // -1 if x=0&z=0, else 0
    // in if x=0&z=0, else out 
    for (int i=0; i<FP2_LIMBS; i++){
        ax[i][0] = jac1[i][0];
        ax[i][1] = jac1[i][2];
        ax[i][2] = jac2[i][0];
        ax[i][3] = jac2[i][2];
    }
    uint32x4_t flag = vandq_u32(theta_point_is_zero(ax), theta_point_is_zero(ax+FP_LIMBS)); //flag = [x, z, x, z]
    az[0] = vdupq_n_u32(flag[0]&flag[1]);
    az[1] = vdupq_n_u32(flag[2]&flag[3]);
    for (int i=0; i<FP_LIMBS; i++){
        //out [x_re, x_im, z_re, z_im]
        //in  [x_re, x_im, z_re, z_im]
        ax[i][0]            = a0[i][0];
        ax[i+FP_LIMBS][0]   = jac1[i][0];
        ax[i][1]            = a0[i+FP_LIMBS][0];
        ax[i+FP_LIMBS][1]   = jac1[i+FP_LIMBS][0];
        ax[i][2]            = a1[i][2];
        ax[i+FP_LIMBS][2]   = jac1[i][2];
        ax[i][3]            = a1[i+FP_LIMBS][2];
        ax[i+FP_LIMBS][3]   = jac1[i+FP_LIMBS][2];
    }
    fp2_select_vec(a2, ax, az[0]);
    for (int i=0; i<FP_LIMBS; i++){
        //x y z
        out1[i][0]          = a2[i][0];
        out1[i+FP_LIMBS][0] = a2[i][1];
        out1[i][1]          = a1[i][0];
        out1[i+FP_LIMBS][1] = a1[i+FP_LIMBS][0];
        out1[i][2]          = a2[i][2];
        out1[i+FP_LIMBS][2] = a2[i][3];
    }

    for (int i=0; i<FP_LIMBS; i++){
        //out [x_re, x_im, z_re, z_im]
        //in  [x_re, x_im, z_re, z_im]
        ax[i][0]          = a0[i][1];
        ax[i+FP_LIMBS][0] = jac2[i][0];
        ax[i][1]          = a0[i+FP_LIMBS][1];
        ax[i+FP_LIMBS][1] = jac2[i+FP_LIMBS][0];
        ax[i][2]          = a1[i][3];
        ax[i+FP_LIMBS][2] = jac2[i][2];
        ax[i][3]          = a1[i+FP_LIMBS][3];
        ax[i+FP_LIMBS][3] = jac2[i+FP_LIMBS][2];
    }
    fp2_select_vec(a2, ax, az[1]);
    for (int i=0; i<FP_LIMBS; i++){
        //x y z
        out2[i][0]          = a2[i][0];
        out2[i+FP_LIMBS][0] = a2[i][1];
        out2[i][1]          = a1[i][1];
        out2[i+FP_LIMBS][1] = a1[i+FP_LIMBS][1];
        out2[i][2]          = a2[i][2];
        out2[i+FP_LIMBS][2] = a2[i][3];
    }
}


void DBLW_vec(uint32x4_t *out)
{ // Cost of 3M + 5S.
  // Doubling on a Weierstrass curve, representation in modified Jacobian coordinates
  // (X:Y:Z:T=a*Z^4) corresponding to (X/Z^2,Y/Z^3), where a is the curve coefficient.
  // Formula from https://hyperelliptic.org/EFD/g1p/auto-shortw-modified.html

    uint32x4_t ax[FP2_LIMBS], az[FP2_LIMBS], a1[FP2_LIMBS], a2[FP2_LIMBS], a3[FP2_LIMBS], TR[FP2_LIMBS], tmp[FP2_LIMBS];

    // uint32_t flag = fp2_is_zero(&in->P1.x) & fp2_is_zero(&in->P1.z);
    // uint32_t flag2 = fp2_is_zero(&in->P2.x) & fp2_is_zero(&in->P2.z);
    for (int i=0; i<FP_LIMBS; i++){
        tmp[i][0] = out[i][0];
        tmp[i][1] = out[i+FP2_LIMBS][1];
        tmp[i][2] = out[i][2];
        tmp[i][3] = out[i+FP2_LIMBS][3];
        tmp[i+FP_LIMBS][0] = out[i+FP_LIMBS][0];
        tmp[i+FP_LIMBS][1] = out[i+FP2_LIMBS+FP_LIMBS][1];
        tmp[i+FP_LIMBS][2] = out[i+FP_LIMBS][2];
        tmp[i+FP_LIMBS][3] = out[i+FP2_LIMBS+FP_LIMBS][3];
    }
    uint32x4_t flag = vandq_u32(theta_point_is_zero(tmp), theta_point_is_zero(tmp+FP_LIMBS));


    // fp2_sqr(&xx, &out->P1.x);               //xx = x^2   [4y+3] 
    // fp2_sqr(&c, &out->P1.y);
    // fp2_sqr(&xx2, &out->P2.x);              //xx = x^2   [4y+3]
    // fp2_sqr(&c2, &out->P2.y);
    fp2_sqr_batched(a1, out);

    
    // fp2_add(&m, &xx, &xx);                   //m = 2x^2   [4y+3]
    // fp2_add(&c, &c, &c);                     //c = 2y^2   [2y+1]
    // fp2_add(&m2, &xx2, &xx2);                //m = 2x^2   [4y+3]
    // fp2_add(&c2, &c2, &c2);                  //c = 2y^2   [2y+1]
    fp2_add_batched(a2, a1, a1);
    fp2_bactched_reduction(a2);


    // fp2_add(&m, &m, &t[0]);                  //xx = 2x^2+t [4y+3]
    // fp2_add(&d, &c, &c);                     //d = 4y^2    [2y+1]
    // fp2_add(&m2, &m2, &t[1]);                //xx = 2x^2+t [4y+3]
    // fp2_add(&d2, &c2, &c2);                  //d = 4y^2    [2y+1]
    /*It must requires to do reduction since 4y^2 is 4 times 29-bit*/
    for (int i=0; i<FP2_LIMBS; i++){
        a3[i][0] = (out+FP2_LIMBS)[i][0];
        a3[i][2] = (out+FP2_LIMBS)[i][2];
        a3[i][1] = a2[i][1];
        a3[i][3] = a2[i][3];
    }
    fp2_add_batched(a3, a2, a3);
    fp2_bactched_reduction(a3);

    // fp2_mul(&d, &d, &out->P1.x);                    //d = 4y^2x  [4y+3]
    // fp2_mul(&out->P1.z, &out->P1.y, &out->P1.z);    //r = yz     [2y+1]
    // fp2_mul(&d2, &d2, &out->P2.x);                  //d = 4y^2x  [4y+3]
    // fp2_mul(&out->P2.z, &out->P2.y, &out->P2.z);    //r = yz     [2y+1]
    for (int i=0; i<FP2_LIMBS; i++){
        a2[i][0] = a3[i][0];
        a2[i][2] = a3[i][2];
        a3[i][0] = a3[i][1];
        a3[i][2] = a3[i][3];
        a3[i][1] = (out+FP2_LIMBS)[i][1];
        a3[i][3] = (out+FP2_LIMBS)[i][3];
    }
    fp2_mul_batched(ax, a3, out);
    

    // fp2_add(&m, &m, &xx);                           //m = 3x^2+t [4y+3]
    // fp2_add(&out->P1.z, &out->P1.z, &out->P1.z);    //outz = 2yz [2y+1]
    // fp2_add(&m2, &m2, &xx2);                        //m = 3x^2+t [4y+3]
    // fp2_add(&out->P2.z, &out->P2.z, &out->P2.z);    //outz = 2yz [2y+1]
    for (int i=0; i<FP2_LIMBS; i++){
        a1[i][1] = a3[i][1] = ax[i][1];                //a1 = [x^2, yz, x^2, yz]
        a1[i][3] = a3[i][3] = ax[i][3];                //a3 = [2x^2+t, yz, 2x^2+t, yz]
        a3[i][0] = a2[i][0];
        a3[i][2] = a2[i][2];
    }
    fp2_add_batched(TR, a1, a3);
    fp2_bactched_reduction(TR);


    // fp2_sqr(&out->P1.x, &m);                    //outx = m^2 [8y+7]
    // fp2_sqr(&cc, &c);                           //cc = 4y^4  [4y+3]
    // fp2_sqr(&out->P2.x, &m2);                   //outx = m^2 [8y+7]
    // fp2_sqr(&cc2, &c2);                         //cc = 4y^4  [4y+3]
    for (int i=0; i<FP2_LIMBS; i++){
        a2[i][0] = TR[i][0];
        a2[i][2] = TR[i][2];
    }
    fp2_sqr_batched(a2, a2);


    // fp2_add(&s, &d, &d);                    //s = 8y^2x     [4y+3]
    // fp2_add(&out->P1.y, &cc, &cc);          //outy = 8y^4   [4y+3]
    // fp2_add(&s2, &d2, &d2);                 //s = 8y^2x     [4y+3]
    // fp2_add(&out->P2.y, &cc2, &cc2);        //outy = 8y^4   [4y+3]
    for (int i=0; i<FP2_LIMBS; i++){
        ax[i][1] = a2[i][1];
        ax[i][3] = a2[i][3];
        a2[i][1] = ax[i][0];
        a2[i][3] = ax[i][2];
    }
    fp2_add_batched(ax, ax, ax);
    fp2_bactched_reduction(ax);

    

    // fp2_sub(&s, &s, &out->P1.x);                //s = 8y^2x - m^2    [4y+3]
    // fp2_sub(&out->P1.x, &out->P1.x, &d);        //outx = m^2 - 4y^2x 
    // fp2_sub(&s2, &s2, &out->P2.x);              //s = 8y^2x - m^2 [4y+3]
    // fp2_sub(&out->P2.x, &out->P2.x, &d2);       //outx = m^2 - 4y^2x
    for (int i=0; i<FP2_LIMBS; i++){
        a3[i][0] = ax[i][0];
        a3[i][2] = ax[i][2];
        a3[i][1] = a2[i][0];
        a3[i][3] = a2[i][2];
    }
    fp2_sub_batched(a3, a3, a2);


    // fp2_add(&t[0], &t[0], &t[0]);               //t = 2t    [4y+3]
    // fp2_add(&s, &s, &d);                        //s = 12y^2x - m^2 [4y+3]
    // fp2_add(&t[1], &t[1], &t[1]);               //t = 2t    [4y+3]
    // fp2_add(&s2, &s2, &d2);                     //s = 12y^2x - m^2 [4y+3]
    for (int i=0; i<FP2_LIMBS; i++){
        a2[i][0]  = (out+FP2_LIMBS)[i][0];
        a2[i][2]  = (out+FP2_LIMBS)[i][2];
        tmp[i][0] = (out+FP2_LIMBS)[i][0];
        tmp[i][1] = a3[i][0];
        tmp[i][2] = (out+FP2_LIMBS)[i][2];
        tmp[i][3] = a3[i][2];
    }
    fp2_add_batched(az, a2, tmp);
    fp2_bactched_reduction(az);


    // fp2_mul(&t[0], &t[0], &out->P1.y);          //t = 2t(8y^4) [8y+7]
    // fp2_mul(&m, &m, &s);                        //m = m(12y^2x - m^2)   [8y+7]
    // fp2_mul(&t[1], &t[1], &out->P2.y);          //t = 2t(8y^4) [8y+7]
    // fp2_mul(&m2, &m2, &s2);                     //m = m(12y^2x - m^2)   [8y+7]
    for (int i=0; i<FP2_LIMBS; i++){
        a1[i][1] = TR[i][0];
        a1[i][3] = TR[i][2];
        a1[i][0] = ax[i][1];
        a1[i][2] = ax[i][3];
    }
    fp2_mul_batched(az, a1, az);
    for (int i=0; i<FP2_LIMBS; i++){                  //TR = [t, z, t, z]
        TR[i][0] = az[i][0];
        TR[i][2] = az[i][2];
    }
    

    // fp2_sub(&out->P1.x, &out->P1.x, &d);    //outx = (m^2 - 4y^2x) - 4y^2x [8y+7]
    // fp2_sub(&out->P1.y, &m, &out->P1.y);    //outy = m(12y^2x - m^2) - (8y^4) [8y+7]
    // fp2_sub(&out->P2.x, &out->P2.x, &d2);   //outx = (m^2 - 4y^2x) - 4y^2x [8y+7]
    // fp2_sub(&out->P2.y, &m2, &out->P2.y);   //outy = m(12y^2x - m^2) - (8y^4) [8y+7]
    for (int i=0; i<FP2_LIMBS; i++){
        az[i][0] = a3[i][1];
        az[i][2] = a3[i][3];
        ax[i][0] = a2[i][1];
        ax[i][2] = a2[i][3];
    }
    fp2_sub_batched(tmp, az, ax);           //tmp = [x, y, x, y]
    fp2_bactched_reduction(tmp);


    // fp2_select(&Q->x, &Q->x, &P->x, -flag);
    // fp2_select(&Q->z, &Q->z, &P->z, -flag);
    // flag = [x, z, x, z]
    // out1 = [x, y, x, y]; out2 = [t, z, t, z]
    az[0] = vdupq_n_u32(flag[0]&flag[1]);
    az[1] = vdupq_n_u32(flag[2]&flag[3]);
    for (int i=0; i<FP_LIMBS; i++){
        //out(Q) [x_re, x_im, z_re, z_im]
        //in(P)  [x_re, x_im, z_re, z_im]
        ax[i][0]            = tmp[i][0];
        ax[i+FP_LIMBS][0]   = out[i][0];
        ax[i][1]            = tmp[i+FP_LIMBS][0];
        ax[i+FP_LIMBS][1]   = out[i+FP_LIMBS][0];
        ax[i][2]            = TR[i][1];
        ax[i+FP_LIMBS][2]   = out[i+FP2_LIMBS][1];
        ax[i][3]            = TR[i+FP_LIMBS][1];
        ax[i+FP_LIMBS][3]   = out[i+FP2_LIMBS+FP_LIMBS][1];
    }
    fp2_select_vec(a2, ax, az[0]);
    for (int i=0; i<FP_LIMBS; i++){
        out[i][0]           = a2[i][0];
        out[i+FP_LIMBS][0]  = a2[i][1];
        out[i][1]           = tmp[i][1];
        out[i+FP_LIMBS][1]  = tmp[i+FP_LIMBS][1];

        out[i+FP2_LIMBS][0]          = TR[i][0];
        out[i+FP2_LIMBS+FP_LIMBS][0] = TR[i+FP_LIMBS][0];
        out[i+FP2_LIMBS][1]          = a2[i][2];
        out[i+FP2_LIMBS+FP_LIMBS][1] = a2[i][3];
    }

    // fp2_select(&Q->x, &Q->x, &P->x, -flag);
    // fp2_select(&Q->z, &Q->z, &P->z, -flag);
    for (int i=0; i<FP_LIMBS; i++){
        //out(Q) [x_re, x_im, z_re, z_im]
        //in(P)  [x_re, x_im, z_re, z_im]
        ax[i][0]            = tmp[i][2];
        ax[i+FP_LIMBS][0]   = out[i][2];
        ax[i][1]            = tmp[i+FP_LIMBS][2];
        ax[i+FP_LIMBS][1]   = out[i+FP_LIMBS][2];
        ax[i][2]            = TR[i][3];
        ax[i+FP_LIMBS][2]   = out[i+FP2_LIMBS][3];
        ax[i][3]            = TR[i+FP_LIMBS][3];
        ax[i+FP_LIMBS][3]   = out[i+FP2_LIMBS+FP_LIMBS][3];
    }
    fp2_select_vec(a2, ax, az[1]);
    for (int i=0; i<FP_LIMBS; i++){
        out[i][2]          = a2[i][0];
        out[i+FP_LIMBS][2] = a2[i][1];
        out[i][3]          = tmp[i][3];
        out[i+FP_LIMBS][3] = tmp[i+FP_LIMBS][3];

        out[i+FP2_LIMBS][2]          = TR[i][2];
        out[i+FP2_LIMBS+FP_LIMBS][2] = TR[i+FP_LIMBS][2];
        out[i+FP2_LIMBS][3]          = a2[i][2];
        out[i+FP2_LIMBS+FP_LIMBS][3] = a2[i][3];
    }
}


void double_couple_jac_point_vec(uint32x4_t *out1, uint32x4_t *out2,
                        const uint32x4_t *jac1,
                        const uint32x4_t *jac2,
                        const uint32x4_t *E1, const uint32x4_t *E2)
{
    DBL_vec(out1, out2, jac1, jac2, E1, E2);
}


void
double_couple_jac_point_iter(theta_couple_jac_point_t *out,
                             unsigned n,
                             const theta_couple_jac_point_t *in,
                             const theta_couple_curve_t *E1E2)
{
    //uint64_t time;
    if (n == 0) {
        *out = *in;
    } else if (n == 1) {
        double_couple_jac_point(out, in, E1E2);
    } else {
        fp2_t a1, a2, t1, t2;

        //time = rdtsc();
        jac_to_ws(&out->P1, &t1, &a1, &in->P1, &E1E2->E1);
        jac_to_ws(&out->P2, &t2, &a2, &in->P2, &E1E2->E2);
        //printf("- jac_to_ws: %lu\n", rdtsc()-time);

        //time = rdtsc();
        DBLW(&out->P1, &t1, &out->P1, &t1);
        DBLW(&out->P2, &t2, &out->P2, &t2);
        for (unsigned i = 0; i < n - 1; i++) {
            DBLW(&out->P1, &t1, &out->P1, &t1);
            DBLW(&out->P2, &t2, &out->P2, &t2);
        }
        //printf("- DBLW: %lu\n", rdtsc()-time);

        //time = rdtsc();
        jac_from_ws(&out->P1, &out->P1, &a1, &E1E2->E1);
        jac_from_ws(&out->P2, &out->P2, &a2, &E1E2->E2);
        //printf("- jac_from_ws: %lu\n", rdtsc()-time);
    }
}


void
jac_to_ws_vec(uint32x4_t *out1, uint32x4_t *out2, uint32x4_t *ta, 
                        const uint32x4_t *jac1, const uint32x4_t *jac2,
                        const uint32x4_t *E1, const uint32x4_t *E2)
{
    // Do two-set jac_to_ws
    // Cost of 3M + 2S when A != 0.
    /* a = 1 - A^2/3, U = X + (A*Z^2)/3, V = Y, W = Z, T = a*Z^4*/

    // fp_div3(&(ao3->re), &(curve1->A.re));    // A/3
    // fp_div3(&(ao3->im), &(curve1->A.im));
    // fp_div3(&(ao3->re), &(curve2->A.re));
    // fp_div3(&(ao3->im), &(curve2->A.im));
    uint32_t inv3[] = THREE_INV_VEC;
    uint32x4_t ao3[FP2_LIMBS] = {0}, tmp[FP2_LIMBS], ax[FP2_LIMBS];
    for (int i=0; i<FP_LIMBS; i++){
        ao3[i][0] = E1[i][0];
        ao3[i][1] = E1[i+FP_LIMBS][0];
        ao3[i][2] = E2[i][0];
        ao3[i][3] = E2[i+FP_LIMBS][0];
        tmp[i] = vdupq_n_u32(inv3[i]);
    }
    fp_mul_batched((uint32x2_t*)ao3, ao3, tmp);    //ao3 = [A/3, A/3, 0, 0]
    for (int i=0; i<FP_LIMBS; i++){
        ao3[i+FP_LIMBS][0] = ao3[i][1];
        ao3[i][1]   = ao3[i][2];
        ao3[i+FP_LIMBS][1] = ao3[i][3];
        ao3[i][2] = ao3[i+FP_LIMBS][2] = ao3[i][3] = ao3[i+FP_LIMBS][3] = 0;
    }

    // fp2_sqr(t, &P->z);                      // t  = z^2
    // fp2_sqr(t, &P->z);                      // t  = z^2
    uint32x4_t t[FP2_LIMBS];
    for (int i=0; i<FP2_LIMBS; i++){
        t[i][0] = jac1[i][2];
        t[i][1] = jac2[i][2];
        t[i][2] = t[i][3] = 0;
    }
    fp2_sqr_batched(t, t);           //t = [z^2 , z^2, 0, 0]

    // fp2_mul(&Q->x, ao3, t);                 // Qx = Az^2/3
    // fp2_mul(&Q->x, ao3, t);                 // Qx = Az^2/3
    fp2_mul_batched(tmp, ao3, t);    //tmp = [Az^2/3, Az^2/3, 0, 0]

    // fp2_add(&Q->x, &Q->x, &P->x);           // Qx = Px + Az^2/3
    // fp2_add(&Q->x, &Q->x, &P->x);           // Qx = Px + Az^2/3
    for (int i=0; i<FP2_LIMBS; i++){
        ax[i][0] = jac1[i][0];
        ax[i][1] = jac2[i][0];
        ax[i][2] = ax[i][3] = 0;
    }
    fp2_add_batched(tmp, tmp, ax);   //tmp = [x+Az^2/3, x+Az^2/3, 0, 0]
    fp2_bactched_reduction(tmp);

    // fp2_sqr(t, t);                          // t  = z^4
    // fp2_sqr(t, t);                          // t  = z^4
    // fp2_mul(&a, ao3, &(curve->A));          // a  = A^2/3
    // fp2_mul(&a, ao3, &(curve->A));          // a  = A^2/3
    for (int i=0; i<FP2_LIMBS; i++){
        t[i][2] = ao3[i][0];
        t[i][3] = ao3[i][1];

        ax[i][0] = t[i][0];
        ax[i][1] = t[i][1];
        ax[i][2] = E1[i][0];
        ax[i][3] = E2[i][0];
    }
    //t = [t, t, ao3, ao3], ax = [t, t, &(curve1->A), &(curve2->A)]
    fp2_mul_batched(ta, t, ax);

    // fp_sub(&(a.re), &one, &(a.re));         // a_re  = 1 - A^2/3
    // fp_neg(&(a.im), &(a.im));               // a_im  = 0 - A^2/3
    // fp_sub(&(a.re), &one, &(a.re));         // a_re  = 1 - A^2/3
    // fp_neg(&(a.im), &(a.im));               // a_im  = 0 - A^2/3
    //q0 = q, q1 = (2^{PER_LIMB*FP_LIMBS}+q)%q
    uint32_t q0[FP_LIMBS], q1[FP_LIMBS];
    for (int i=0; i<(FP_LIMBS-1); i++) q0[i] = Q1_VALUE;
    q0[FP_LIMBS-1] = Q1_VALUE_SIGNIFICANT;
    initial_q_value_1(q1);
    
    for (int i=0; i<FP_LIMBS; i++){
        t[i][0] = q1[i];
        t[i][1] = q0[i];
        t[i][2] = q1[i];
        t[i][3] = q0[i];

        ax[i][0] = ta[i][2];
        ax[i][1] = ta[i+FP_LIMBS][2];
        ax[i][2] = ta[i][3];
        ax[i][3] = ta[i+FP_LIMBS][3];
    }
    fp_sub_batched(t, t, ax);                 // t = [1 - A^2/3, 0 - A^2/3, 1 - A^2/3, 0 - A^2/3]
    fp_bactched_reduction(t);
    
    // fp2_mul(t, t, &a);                      // t_re = (t_re*a_re) - (t_im*a_im); t_im = (t_im*a_re) - (t_re*a_im)
    // fp2_mul(t, t, &a);                      // t_re = (t_re*a_re) - (t_im*a_im); t_im = (t_im*a_re) - (t_re*a_im)
    for (int i=0; i<FP_LIMBS; i++){
        ax[i][0]   = t[i][0];
        ax[i+FP_LIMBS][0] = t[i][1];
        ax[i][1]   = t[i][2];
        ax[i+FP_LIMBS][1] = t[i][3];
        ax[i][2] = ax[i+FP_LIMBS][2] = ax[i][3] = ax[i+FP_LIMBS][3] = 0;
    }
    fp2_mul_batched(ta, ta, ax);     //ta = [t1, t2, a1, a2]
    for (int i=0; i<FP2_LIMBS; i++){
        ta[i][2] = ao3[i][0];
        ta[i][3] = ao3[i][1];
    }

    // fp2_copy(&Q->y, &P->y);
    // fp2_copy(&Q->z, &P->z);
    for (int i=0; i<FP2_LIMBS; i++){
        out1[i][0] = tmp[i][0];
        out1[i][1] = jac1[i][1];
        out1[i][2] = jac1[i][2];

        out2[i][0] = tmp[i][1];
        out2[i][1] = jac2[i][1];
        out2[i][2] = jac2[i][2];
    }
}


void
jac_from_ws_vec(uint32x4_t *out1, uint32x4_t *out2, const uint32x4_t *in1, const uint32x4_t *in2, const uint32x4_t *ta, const uint32x4_t *E1, const uint32x4_t *E2)
{
    // Cost of 1M + 1S when A != 0.
    /* X = U - (A*W^2)/3, Y = V, Z = W. */
    
    // fp2_sqr(&t, &P->z);
    // fp2_sqr(&t, &P->z);
    uint32x4_t tmp[FP2_LIMBS], ao3[FP2_LIMBS];
    for (int i=0; i<FP2_LIMBS; i++){
        tmp[i][0] = in1[i][2];
        tmp[i][1] = in2[i][2];
        tmp[i][2] = tmp[i][3] = 0;

        ao3[i][0] = ta[i][2];
        ao3[i][1] = ta[i][3];
        ao3[i][2] = ao3[i][3] = 0;
    }
    fp2_sqr_batched(tmp, tmp);

    // fp2_mul(&t, &t, ao3);
    // fp2_mul(&t, &t, ao3);
    fp2_mul_batched(ao3, tmp, ao3);

    // fp2_sub(&Q->x, &P->x, &t);
    // fp2_sub(&Q->x, &P->x, &t);
    for (int i=0; i<FP2_LIMBS; i++){
        tmp[i][0] = in1[i][0];
        tmp[i][1] = in2[i][0];
        tmp[i][2] = tmp[i][3] = 0;
    }
    fp2_sub_batched(tmp, tmp, ao3);
    fp2_bactched_reduction(tmp);
    
    // fp2_copy(&Q->y, &P->y);
    // fp2_copy(&Q->z, &P->z);
    // fp2_copy(&Q->y, &P->y);
    // fp2_copy(&Q->z, &P->z);
    for (int i=0; i<FP2_LIMBS; i++){
        out1[i][0] = tmp[i][0];
        out1[i][1] = in1[i][1];
        out1[i][2] = in1[i][2];

        out2[i][0] = tmp[i][1];
        out2[i][1] = in2[i][1];
        out2[i][2] = in2[i][2];
    }

}


//jac1 = [x, y, z, non-def]; jac2 = [x, y, z, non-def]
//E1 = [A, C, A24.x, A24.z]; E2 = [A, C, A24.x, A24.z]
void double_couple_jac_point_iter_vec(uint32x4_t *out1, uint32x4_t *out2,
                        unsigned n,
                        const uint32x4_t *jac1, const uint32x4_t *jac2,
                        const uint32x4_t *E1, const uint32x4_t *E2)
{
    if (n == 0) {
        memcpy(out1, jac1, FP2_LIMBS*sizeof(uint32x4_t));
        memcpy(out2, jac2, FP2_LIMBS*sizeof(uint32x4_t));
    } else if (n == 1) {
        double_couple_jac_point_vec(out1, out2, jac1, jac2, E1, E2);
    } else {
        // fp2_t a1, a2, t1, t2;
        // jac_to_ws(&out->P1, &t1, &a1, &in->P1, &E1E2->E1);
        // jac_to_ws(&out->P2, &t2, &a2, &in->P2, &E1E2->E2);
        uint32x4_t ta[FP2_LIMBS]; //[t1, t2, a1, a2]
        jac_to_ws_vec(out1, out2, ta, jac1, jac2, E1, E2);

        // DBLW(&out->P2, &t2, &out->P2, &t2);
        // for (unsigned i = 0; i < n - 1; i++){
        //     DBLW(&out->P1, &t1, &out->P1, &t1);
        //     DBLW(&out->P2, &t2, &out->P2, &t2);
        // }
        // [x, y, x, y]; [t, z, t, z]
        uint32x4_t tmp[FP2_LIMBS*2];
        for (int i=0; i<FP2_LIMBS; i++){
            tmp[i][0] = out1[i][0];
            tmp[i][1] = out1[i][1];
            tmp[i][2] = out2[i][0];
            tmp[i][3] = out2[i][1];

            tmp[i+FP2_LIMBS][0] = ta[i][0];
            tmp[i+FP2_LIMBS][1] = out1[i][2];
            tmp[i+FP2_LIMBS][2] = ta[i][1];
            tmp[i+FP2_LIMBS][3] = out2[i][2];
        }
        DBLW_vec(tmp);
        for (unsigned i = 0; i < n - 1; i++) {
            DBLW_vec(tmp);
        }
        //correct position: x y z
        for (int i=0; i<FP2_LIMBS; i++){
            out1[i][0] = tmp[i][0];
            out1[i][1] = tmp[i][1];
            out1[i][2] = tmp[i+FP2_LIMBS][1];

            out2[i][0] = tmp[i][2];
            out2[i][1] = tmp[i][3];
            out2[i][2] = tmp[i+FP2_LIMBS][3];
        }

        // jac_from_ws(&out->P1, &out->P1, &a1, &E1E2->E1);
        // jac_from_ws(&out->P2, &out->P2, &a2, &E1E2->E2);
        jac_from_ws_vec(out1, out2, out1, out2, ta, E1, E2);
    }
}

void
couple_jac_to_xz(theta_couple_point_t *P, const theta_couple_jac_point_t *xyP)
{
    jac_to_xz(&P->P1, &xyP->P1);
    jac_to_xz(&P->P2, &xyP->P2);
}

void
copy_bases_to_kernel(theta_kernel_couple_points_t *ker, const ec_basis_t *B1, const ec_basis_t *B2)
{
    // Copy the basis on E1 to (P, _) on T1, T2 and T1 - T2
    copy_point(&ker->T1.P1, &B1->P);
    copy_point(&ker->T2.P1, &B1->Q);
    copy_point(&ker->T1m2.P1, &B1->PmQ);

    // Copy the basis on E2 to (_, P) on T1, T2 and T1 - T2
    copy_point(&ker->T1.P2, &B2->P);
    copy_point(&ker->T2.P2, &B2->Q);
    copy_point(&ker->T1m2.P2, &B2->PmQ);
}

int xDBLMUL_vec(ec_point_t *S,
        const ec_point_t *P,
        const digit_t *k,
        const ec_point_t *Q,
        const digit_t *l,
        const ec_point_t *PQ,
        const int kbits,
        const ec_curve_t *curve)                    
{
    // The Montgomery biladder
    // Input:  projective Montgomery points P=(XP:ZP) and Q=(XQ:ZQ) such that xP=XP/ZP and xQ=XQ/ZQ, scalars k and l of
    //         bitlength kbits, the difference PQ=P-Q=(XPQ:ZPQ), and the Montgomery curve constants (A:C).
    // Output: projective Montgomery point S <- k*P + l*Q = (XS:ZS) such that x(k*P + l*Q)=XS/ZS.
    int i, A_is_zero;
    digit_t evens, mevens, bitk0, bitl0, maskk, maskl, temp, bs1_ip1, bs2_ip1, bs1_i, bs2_i, h;
    digit_t sigma[2] = { 0 }, pre_sigma = 0;
    digit_t k_t[NWORDS_ORDER], l_t[NWORDS_ORDER], one[NWORDS_ORDER] = { 0 }, r[2 * BITS] = { 0 };
    ec_point_t R[3] = {0};

    // differential additions formulas are invalid in this case
    if (ec_has_zero_coordinate(P) | ec_has_zero_coordinate(Q) | ec_has_zero_coordinate(PQ))
        return 0;

    // Derive sigma according to parity
    bitk0 = (k[0] & 1);
    bitl0 = (l[0] & 1);
    maskk = 0 - bitk0; // Parity masks: 0 if even, otherwise 1...1
    maskl = 0 - bitl0;
    sigma[0] = (bitk0 ^ 1);
    sigma[1] = (bitl0 ^ 1);
    evens = sigma[0] + sigma[1]; // Count number of even scalars
    mevens = 0 - (evens & 1);    // Mask mevens <- 0 if # even of scalars = 0 or 2, otherwise mevens = 1...1

    // If k and l are both even or both odd, pick sigma = (0,1)
    sigma[0] = (sigma[0] & mevens);
    sigma[1] = (sigma[1] & mevens) | (1 & ~mevens);

    // Convert even scalars to odd
    one[0] = 1;
    mp_sub(k_t, k, one, NWORDS_ORDER);
    mp_sub(l_t, l, one, NWORDS_ORDER);

    select_ct(k_t, k_t, k, maskk, NWORDS_ORDER);
    select_ct(l_t, l_t, l, maskl, NWORDS_ORDER);


    // Scalar recoding
    for (i = 0; i < kbits; i++) {
        // If sigma[0] = 1 swap k_t and l_t
        maskk = 0 - (sigma[0] ^ pre_sigma);
        swap_ct(k_t, l_t, maskk, NWORDS_ORDER);

        if (i == kbits - 1) {
            bs1_ip1 = 0;
            bs2_ip1 = 0;
        } else {
            bs1_ip1 = mp_shiftr(k_t, 1, NWORDS_ORDER);
            bs2_ip1 = mp_shiftr(l_t, 1, NWORDS_ORDER);
        }
        bs1_i = k_t[0] & 1;
        bs2_i = l_t[0] & 1;

        r[2 * i] = bs1_i ^ bs1_ip1;
        r[2 * i + 1] = bs2_i ^ bs2_ip1;

        // Revert sigma if second bit, r_(2i+1), is 1
        pre_sigma = sigma[0];
        maskk = 0 - r[2 * i + 1];
        select_ct(&temp, &sigma[0], &sigma[1], maskk, 1);
        select_ct(&sigma[1], &sigma[1], &sigma[0], maskk, 1);
        sigma[0] = temp;
    }

    // Point initialization
    ec_point_init(&R[0]);
    maskk = 0 - sigma[0];
    select_point(&R[1], P, Q, maskk);
    select_point(&R[2], Q, P, maskk);

    theta_point_t tp;
    fp_t mb, imb;
    initial_montback_array(&mb, &imb);
    uint32x4_t In32[FP2_LIMBS], DIFF1a32[FP_LIMBS], DIFF1b32[FP_LIMBS], DIFF2a32[FP_LIMBS], DIFF2b32[FP_LIMBS];
    tp.x = R[1].x;
    tp.y = R[1].z;
    tp.z = R[2].x;
    tp.t = R[2].z;
    theta_montback(&tp, &mb);
    transpose(In32, tp);
    for (i=0; i<FP_LIMBS; i++){
        DIFF1a32[i][0] = In32[i][0];
        DIFF1a32[i][1] = In32[i+FP_LIMBS][0];
        DIFF1a32[i][2] = In32[i][1];
        DIFF1a32[i][3] = In32[i+FP_LIMBS][1];

        DIFF1b32[i][0] = In32[i][2];
        DIFF1b32[i][1] = In32[i+FP_LIMBS][2];
        DIFF1b32[i][2] = In32[i][3];
        DIFF1b32[i][3] = In32[i+FP_LIMBS][3];
    }

    // Initialize DIFF2a <- P+Q, DIFF2b <- P-Q
    xADD(&R[2], &R[1], &R[2], PQ);

    if (ec_has_zero_coordinate(&R[2]))
        return 0; // non valid formulas

    tp.x = R[2].x;
    tp.y = R[2].z;
    tp.z = PQ->x;
    tp.t = PQ->z;
    theta_montback(&tp, &mb);
    transpose(In32, tp);
    for (i=0; i<FP_LIMBS; i++){
        DIFF2a32[i][0] = In32[i][0];
        DIFF2a32[i][1] = In32[i+FP_LIMBS][0];
        DIFF2a32[i][2] = In32[i][1];
        DIFF2a32[i][3] = In32[i+FP_LIMBS][1];

        DIFF2b32[i][0] = In32[i][2];
        DIFF2b32[i][1] = In32[i+FP_LIMBS][2];
        DIFF2b32[i][2] = In32[i][3];
        DIFF2b32[i][3] = In32[i+FP_LIMBS][3];
    }

    A_is_zero = fp2_is_zero(&curve->A);

    uint32x4_t R01[FP2_LIMBS], R12[FP2_LIMBS], A24[FP2_LIMBS], Arith1[FP2_LIMBS], Arith2[FP2_LIMBS];
    tp.x = R[0].x;
    tp.y = R[0].z;
    tp.z = R[1].x;
    tp.t = R[1].z;
    theta_montback(&tp, &mb);
    transpose(In32, tp);
    for (i=0; i<FP_LIMBS; i++){
        R01[i][0] = In32[i][0];
        R01[i][1] = In32[i+FP_LIMBS][0];
        R01[i][2] = In32[i][1];
        R01[i][3] = In32[i+FP_LIMBS][1];

        R01[i+FP_LIMBS][0] = In32[i][2];
        R01[i+FP_LIMBS][1] = In32[i+FP_LIMBS][2];
        R01[i+FP_LIMBS][2] = In32[i][3];
        R01[i+FP_LIMBS][3] = In32[i+FP_LIMBS][3];
    }

    tp.x = R[1].x;
    tp.y = R[1].z;
    tp.z = R[2].x;
    tp.t = R[2].z;
    theta_montback(&tp, &mb);
    transpose(In32, tp);
    for (i=0; i<FP_LIMBS; i++){
        R12[i][0] = In32[i][0];
        R12[i][1] = In32[i+FP_LIMBS][0];
        R12[i][2] = In32[i][1];
        R12[i][3] = In32[i+FP_LIMBS][1];

        R12[i+FP_LIMBS][0] = In32[i][2];
        R12[i+FP_LIMBS][1] = In32[i+FP_LIMBS][2];
        R12[i+FP_LIMBS][2] = In32[i][3];
        R12[i+FP_LIMBS][3] = In32[i+FP_LIMBS][3];
    }

    for (int i=0; i<NWORDS_FIELD; i++){
        tp.x.re[i] = curve->A24.x.re[i];
        tp.y.re[i] = curve->A24.x.im[i];
        tp.z.re[i] = curve->A24.z.re[i];
        tp.t.re[i] = curve->A24.z.im[i];
    }
    theta_montback(&tp, &mb);
    transpose(A24, tp);
    // Main loop
    for (i = kbits - 1; i >= 0; i--) {
        h = r[2 * i] + r[2 * i + 1]; // in {0, 1, 2}
        maskk = 0 - (h & 1);
        fp2_select_vec(Arith1, R01, vdupq_n_u32(maskk));

        maskk = 0 - (h >> 1);
        for (int z=0; z<FP_LIMBS; z++){
            Arith1[z+FP_LIMBS] = R12[z+FP_LIMBS];
        }
        fp2_select_vec(Arith1, Arith1, vdupq_n_u32(maskk));

        if (A_is_zero) {
            xDBL_E0_vec(Arith1, Arith1);
        } else {
            //assert(fp2_is_one(&curve->A24.z));
            xDBL_A24_vec(Arith1, Arith1, A24, true);
        }

        maskk = 0 - r[2 * i + 1]; // in {0, 1}
        fp2_select_vec(Arith1+FP_LIMBS, R01, vdupq_n_u32(maskk));
        fp2_select_vec(Arith2, R12, vdupq_n_u32(maskk));

        cswap_points_vec(DIFF1a32, DIFF1b32, maskk);
        
        //Before: Arith1 <- T0T1; Arith2 <- T2;
        for (int z=0; z<FP_LIMBS; z++){
            R12[z] = R01[z];
            //R0 update
            R01[z] = Arith1[z];
            Arith1[z] = Arith2[z];
        }
        xADD2_vec(Arith1+FP_LIMBS, Arith1, DIFF1a32, Arith2, R12, DIFF2a32);

        // If hw (mod 2) = 1 then swap DIFF2a and DIFF2b
        maskk = 0 - (h & 1);
        cswap_points_vec(DIFF2a32, DIFF2b32, maskk);

        // R <- T
        for (int z=0; z<FP_LIMBS; z++){
            //R1, R2 update
            R01[z+FP_LIMBS] =  R12[z] = Arith1[z+FP_LIMBS];
            R12[z+FP_LIMBS] = Arith2[z];
        }
    }

    itranspose(&tp, R01);
    theta_montback(&tp, &imb);
    for (int z=0; z<NWORDS_FIELD; z++){
        R[0].x.re[z] = tp.x.re[z];
        R[0].x.im[z] = tp.y.re[z];
        R[0].z.re[z] = tp.z.re[z];
        R[0].z.im[z] = tp.t.re[z];

        R[1].x.re[z] = tp.x.im[z];
        R[1].x.im[z] = tp.y.im[z];
        R[1].z.re[z] = tp.z.im[z];
        R[1].z.im[z] = tp.t.im[z];
    }

    itranspose(&tp, R12);
    theta_montback(&tp, &imb);
    for (int z=0; z<NWORDS_FIELD; z++){
        R[2].x.re[z] = tp.x.im[z];
        R[2].x.im[z] = tp.y.im[z];
        R[2].z.re[z] = tp.z.im[z];
        R[2].z.im[z] = tp.t.im[z];
    }

    // Output R[evens]
    select_point(S, &R[0], &R[1], mevens);
    maskk = 0 - (bitk0 & bitl0);
    select_point(S, S, &R[2], maskk);
    return 1;
}