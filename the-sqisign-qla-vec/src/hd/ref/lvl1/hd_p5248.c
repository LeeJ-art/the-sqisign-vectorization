#include <stdint.h>
#include <stdio.h>
#include <arm_neon.h>
#include <hd.h>
#include <ec.h>

void initial_montback_array(fp_t* mb, fp_t* imb){
    //R2 = 2^261
    mb[0][0] = 1638;
    mb[0][1] = 0;
    mb[0][2] = 0;
    mb[0][3] = 0;
    mb[0][4] = 35184372088832;
    //2^(255*2-261)
    imb[0][0] = 0;
    imb[0][1] = 0;
    imb[0][2] = 0;
    imb[0][3] = 0;
    imb[0][4] = 35184372088832;
}

void initial_q_value_1(uint32_t* q1){
    uint32_t tmp[FP_LIMBS] = {1638, 0, 0, 0, 0, 0, 0, 0, 131072};
    memmove(q1, tmp, sizeof(uint32_t)*FP_LIMBS);
}

void ec_montback_array(fp_t* mb){
    //2^(261*4 - 255*3)
    mb[0][0] = 429496729;
    mb[0][1] = 0;
    mb[0][2] = 0;
    mb[0][3] = 0;
    mb[0][4] = 52776558133248;
}

void transpose(uint32x4_t *Out, theta_point_t In){
    uint64_t mask = ((uint64_t)1<<29)-1;
    uint32_t in32[4][18] = {0};

    // re-cast
    in32[0][0] = (uint32_t)(mask & In.x.re[0]);
    in32[1][0] = (uint32_t)(mask & In.y.re[0]);
    in32[2][0] = (uint32_t)(mask & In.z.re[0]);
    in32[3][0] = (uint32_t)(mask & In.t.re[0]);

    in32[0][1] = (uint32_t)(mask & ((In.x.re[0]>>29) + (In.x.re[1]<<22)));
    in32[1][1] = (uint32_t)(mask & ((In.y.re[0]>>29) + (In.y.re[1]<<22)));
    in32[2][1] = (uint32_t)(mask & ((In.z.re[0]>>29) + (In.z.re[1]<<22)));
    in32[3][1] = (uint32_t)(mask & ((In.t.re[0]>>29) + (In.t.re[1]<<22)));

    in32[0][2] = (uint32_t)(mask & (In.x.re[1]>>7));
    in32[1][2] = (uint32_t)(mask & (In.y.re[1]>>7));
    in32[2][2] = (uint32_t)(mask & (In.z.re[1]>>7));
    in32[3][2] = (uint32_t)(mask & (In.t.re[1]>>7));

    in32[0][3] = (uint32_t)(mask & ((In.x.re[1]>>36) + (In.x.re[2]<<15)));
    in32[1][3] = (uint32_t)(mask & ((In.y.re[1]>>36) + (In.y.re[2]<<15)));
    in32[2][3] = (uint32_t)(mask & ((In.z.re[1]>>36) + (In.z.re[2]<<15)));
    in32[3][3] = (uint32_t)(mask & ((In.t.re[1]>>36) + (In.t.re[2]<<15)));

    in32[0][4] = (uint32_t)(mask & (In.x.re[2]>>14));
    in32[1][4] = (uint32_t)(mask & (In.y.re[2]>>14));
    in32[2][4] = (uint32_t)(mask & (In.z.re[2]>>14));
    in32[3][4] = (uint32_t)(mask & (In.t.re[2]>>14));

    in32[0][5] = (uint32_t)(mask & ((In.x.re[2]>>43) + (In.x.re[3]<<8)));
    in32[1][5] = (uint32_t)(mask & ((In.y.re[2]>>43) + (In.y.re[3]<<8)));
    in32[2][5] = (uint32_t)(mask & ((In.z.re[2]>>43) + (In.z.re[3]<<8)));
    in32[3][5] = (uint32_t)(mask & ((In.t.re[2]>>43) + (In.t.re[3]<<8)));

    in32[0][6] = (uint32_t)(mask & (In.x.re[3]>>21));
    in32[1][6] = (uint32_t)(mask & (In.y.re[3]>>21));
    in32[2][6] = (uint32_t)(mask & (In.z.re[3]>>21));
    in32[3][6] = (uint32_t)(mask & (In.t.re[3]>>21));

    in32[0][7] = (uint32_t)(mask & ((In.x.re[3]>>50) + (In.x.re[4]<<1)));
    in32[1][7] = (uint32_t)(mask & ((In.y.re[3]>>50) + (In.y.re[4]<<1)));
    in32[2][7] = (uint32_t)(mask & ((In.z.re[3]>>50) + (In.z.re[4]<<1)));
    in32[3][7] = (uint32_t)(mask & ((In.t.re[3]>>50) + (In.t.re[4]<<1)));

    in32[0][8] = (uint32_t)(mask & (In.x.re[4]>>28));
    in32[1][8] = (uint32_t)(mask & (In.y.re[4]>>28));
    in32[2][8] = (uint32_t)(mask & (In.z.re[4]>>28));
    in32[3][8] = (uint32_t)(mask & (In.t.re[4]>>28));

    in32[0][9] = (uint32_t)(mask & (In.x.im[0]));
    in32[1][9] = (uint32_t)(mask & (In.y.im[0]));
    in32[2][9] = (uint32_t)(mask & (In.z.im[0]));
    in32[3][9] = (uint32_t)(mask & (In.t.im[0]));

    in32[0][10] = (uint32_t)(mask & ((In.x.im[0]>>29) + (In.x.im[1]<<22)));
    in32[1][10] = (uint32_t)(mask & ((In.y.im[0]>>29) + (In.y.im[1]<<22)));
    in32[2][10] = (uint32_t)(mask & ((In.z.im[0]>>29) + (In.z.im[1]<<22)));
    in32[3][10] = (uint32_t)(mask & ((In.t.im[0]>>29) + (In.t.im[1]<<22)));

    in32[0][11] = (uint32_t)(mask & (In.x.im[1]>>7));
    in32[1][11] = (uint32_t)(mask & (In.y.im[1]>>7));
    in32[2][11] = (uint32_t)(mask & (In.z.im[1]>>7));
    in32[3][11] = (uint32_t)(mask & (In.t.im[1]>>7));

    in32[0][12] = (uint32_t)(mask & ((In.x.im[1]>>36) + (In.x.im[2]<<15)));
    in32[1][12] = (uint32_t)(mask & ((In.y.im[1]>>36) + (In.y.im[2]<<15)));
    in32[2][12] = (uint32_t)(mask & ((In.z.im[1]>>36) + (In.z.im[2]<<15)));
    in32[3][12] = (uint32_t)(mask & ((In.t.im[1]>>36) + (In.t.im[2]<<15)));

    in32[0][13] = (uint32_t)(mask & (In.x.im[2]>>14));
    in32[1][13] = (uint32_t)(mask & (In.y.im[2]>>14));
    in32[2][13] = (uint32_t)(mask & (In.z.im[2]>>14));
    in32[3][13] = (uint32_t)(mask & (In.t.im[2]>>14));

    in32[0][14] = (uint32_t)(mask & ((In.x.im[2]>>43) + (In.x.im[3]<<8)));
    in32[1][14] = (uint32_t)(mask & ((In.y.im[2]>>43) + (In.y.im[3]<<8)));
    in32[2][14] = (uint32_t)(mask & ((In.z.im[2]>>43) + (In.z.im[3]<<8)));
    in32[3][14] = (uint32_t)(mask & ((In.t.im[2]>>43) + (In.t.im[3]<<8)));

    in32[0][15] = (uint32_t)(mask & (In.x.im[3]>>21));
    in32[1][15] = (uint32_t)(mask & (In.y.im[3]>>21));
    in32[2][15] = (uint32_t)(mask & (In.z.im[3]>>21));
    in32[3][15] = (uint32_t)(mask & (In.t.im[3]>>21));

    in32[0][16] = (uint32_t)(mask & ((In.x.im[3]>>50) + (In.x.im[4]<<1)));
    in32[1][16] = (uint32_t)(mask & ((In.y.im[3]>>50) + (In.y.im[4]<<1)));
    in32[2][16] = (uint32_t)(mask & ((In.z.im[3]>>50) + (In.z.im[4]<<1)));
    in32[3][16] = (uint32_t)(mask & ((In.t.im[3]>>50) + (In.t.im[4]<<1)));

    in32[0][17] = (uint32_t)(mask & (In.x.im[4]>>28));
    in32[1][17] = (uint32_t)(mask & (In.y.im[4]>>28));
    in32[2][17] = (uint32_t)(mask & (In.z.im[4]>>28));
    in32[3][17] = (uint32_t)(mask & (In.t.im[4]>>28));

    // transpose
    for(int i = 0;i<18;i++){
        uint32x4_t tmp = {in32[0][i], in32[1][i], in32[2][i], in32[3][i]};
        Out[i] = tmp;
    }
}

void itranspose(theta_point_t *Out, uint32x4_t *In){
    uint32_t in32[4][18];

    // itranspose
    for(int i = 0;i<4;i++){
        for(int j = 0;j<18;j++){
            in32[i][j] = In[j][i];
        }
    }

    // re-cast
    uint64_t mask = ((uint64_t)1<<51)-1;
    Out->x.re[0] = mask & (in32[0][0] + ((uint64_t)in32[0][1]<<29));
    Out->x.im[0] = mask & (in32[0][9] + ((uint64_t)in32[0][10]<<29));
    Out->y.re[0] = mask & (in32[1][0] + ((uint64_t)in32[1][1]<<29));
    Out->y.im[0] = mask & (in32[1][9] + ((uint64_t)in32[1][10]<<29));
    Out->z.re[0] = mask & (in32[2][0] + ((uint64_t)in32[2][1]<<29));
    Out->z.im[0] = mask & (in32[2][9] + ((uint64_t)in32[2][10]<<29));
    Out->t.re[0] = mask & (in32[3][0] + ((uint64_t)in32[3][1]<<29));
    Out->t.im[0] = mask & (in32[3][9] + ((uint64_t)in32[3][10]<<29));

    Out->x.re[1] = mask & (((uint64_t)in32[0][1]>>22) + ((uint64_t)in32[0][2]<<7) + ((uint64_t)in32[0][3]<<36));
    Out->x.im[1] = mask & (((uint64_t)in32[0][10]>>22) + ((uint64_t)in32[0][11]<<7) + ((uint64_t)in32[0][12]<<36));
    Out->y.re[1] = mask & (((uint64_t)in32[1][1]>>22) + ((uint64_t)in32[1][2]<<7) + ((uint64_t)in32[1][3]<<36));
    Out->y.im[1] = mask & (((uint64_t)in32[1][10]>>22) + ((uint64_t)in32[1][11]<<7) + ((uint64_t)in32[1][12]<<36));
    Out->z.re[1] = mask & (((uint64_t)in32[2][1]>>22) + ((uint64_t)in32[2][2]<<7) + ((uint64_t)in32[2][3]<<36));
    Out->z.im[1] = mask & (((uint64_t)in32[2][10]>>22) + ((uint64_t)in32[2][11]<<7) + ((uint64_t)in32[2][12]<<36));
    Out->t.re[1] = mask & (((uint64_t)in32[3][1]>>22) + ((uint64_t)in32[3][2]<<7) + ((uint64_t)in32[3][3]<<36));
    Out->t.im[1] = mask & (((uint64_t)in32[3][10]>>22) + ((uint64_t)in32[3][11]<<7) + ((uint64_t)in32[3][12]<<36));

    Out->x.re[2] = mask & (((uint64_t)in32[0][3]>>15) + ((uint64_t)in32[0][4]<<14) + ((uint64_t)in32[0][5]<<43));
    Out->x.im[2] = mask & (((uint64_t)in32[0][12]>>15) + ((uint64_t)in32[0][13]<<14) + ((uint64_t)in32[0][14]<<43));
    Out->y.re[2] = mask & (((uint64_t)in32[1][3]>>15) + ((uint64_t)in32[1][4]<<14) + ((uint64_t)in32[1][5]<<43));
    Out->y.im[2] = mask & (((uint64_t)in32[1][12]>>15) + ((uint64_t)in32[1][13]<<14) + ((uint64_t)in32[1][14]<<43));
    Out->z.re[2] = mask & (((uint64_t)in32[2][3]>>15) + ((uint64_t)in32[2][4]<<14) + ((uint64_t)in32[2][5]<<43));
    Out->z.im[2] = mask & (((uint64_t)in32[2][12]>>15) + ((uint64_t)in32[2][13]<<14) + ((uint64_t)in32[2][14]<<43));
    Out->t.re[2] = mask & (((uint64_t)in32[3][3]>>15) + ((uint64_t)in32[3][4]<<14) + ((uint64_t)in32[3][5]<<43));
    Out->t.im[2] = mask & (((uint64_t)in32[3][12]>>15) + ((uint64_t)in32[3][13]<<14) + ((uint64_t)in32[3][14]<<43));

    Out->x.re[3] = mask & (((uint64_t)in32[0][5]>>8) + ((uint64_t)in32[0][6]<<21) + ((uint64_t)in32[0][7]<<50));
    Out->x.im[3] = mask & (((uint64_t)in32[0][14]>>8) + ((uint64_t)in32[0][15]<<21) + ((uint64_t)in32[0][16]<<50));
    Out->y.re[3] = mask & (((uint64_t)in32[1][5]>>8) + ((uint64_t)in32[1][6]<<21) + ((uint64_t)in32[1][7]<<50));
    Out->y.im[3] = mask & (((uint64_t)in32[1][14]>>8) + ((uint64_t)in32[1][15]<<21) + ((uint64_t)in32[1][16]<<50));
    Out->z.re[3] = mask & (((uint64_t)in32[2][5]>>8) + ((uint64_t)in32[2][6]<<21) + ((uint64_t)in32[2][7]<<50));
    Out->z.im[3] = mask & (((uint64_t)in32[2][14]>>8) + ((uint64_t)in32[2][15]<<21) + ((uint64_t)in32[2][16]<<50));
    Out->t.re[3] = mask & (((uint64_t)in32[3][5]>>8) + ((uint64_t)in32[3][6]<<21) + ((uint64_t)in32[3][7]<<50));
    Out->t.im[3] = mask & (((uint64_t)in32[3][14]>>8) + ((uint64_t)in32[3][15]<<21) + ((uint64_t)in32[3][16]<<50));

    Out->x.re[4] = mask & (((uint64_t)in32[0][7]>>1) + ((uint64_t)in32[0][8]<<28));
    Out->x.im[4] = mask & (((uint64_t)in32[0][16]>>1) + ((uint64_t)in32[0][17]<<28));
    Out->y.re[4] = mask & (((uint64_t)in32[1][7]>>1) + ((uint64_t)in32[1][8]<<28));
    Out->y.im[4] = mask & (((uint64_t)in32[1][16]>>1) + ((uint64_t)in32[1][17]<<28));
    Out->z.re[4] = mask & (((uint64_t)in32[2][7]>>1) + ((uint64_t)in32[2][8]<<28));
    Out->z.im[4] = mask & (((uint64_t)in32[2][16]>>1) + ((uint64_t)in32[2][17]<<28));
    Out->t.re[4] = mask & (((uint64_t)in32[3][7]>>1) + ((uint64_t)in32[3][8]<<28));
    Out->t.im[4] = mask & (((uint64_t)in32[3][16]>>1) + ((uint64_t)in32[3][17]<<28));
}