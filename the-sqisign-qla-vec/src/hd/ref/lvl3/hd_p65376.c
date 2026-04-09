#include <stdint.h>
#include <stdio.h>
#include <arm_neon.h>
#include <hd.h>

void initial_montback_array(fp_t* mb, fp_t* imb){
    //R2 = 2^392
    mb[0][0] = 1008;
    mb[0][1] = 0;
    mb[0][2] = 0;
    mb[0][3] = 0;
    mb[0][4] = 0;
    mb[0][5] = 0;
    mb[0][6] = 1125899906842624;
    //2^(385*2-392)
    imb[0][0] = 0;
    imb[0][1] = 0;
    imb[0][2] = 0;
    imb[0][3] = 0;
    imb[0][4] = 0;
    imb[0][5] = 0;
    imb[0][6] = 281474976710656;
}

void initial_q_value_1(uint32_t* q1){
    uint32_t tmp[FP_LIMBS] = {1008, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 65536};
    memmove(q1, tmp, sizeof(uint32_t)*FP_LIMBS);
}

void ec_montback_array(fp_t* mb){
    //2^(392*4 - 385*3)
    mb[0][0] = 2114445438;
    mb[0][1] = 0;
    mb[0][2] = 0;
    mb[0][3] = 0;
    mb[0][4] = 0;
    mb[0][5] = 0;
    mb[0][6] = 140737488355328;
}

void transpose(uint32x4_t *Out, theta_point_t In){
    uint64_t mask = ((uint64_t)1<<PER_LIMB)-1;
    uint32_t in32[4][FP2_LIMBS] = {0};

    // re-cast
    in32[0][0] = (uint32_t)(mask & In.x.re[0]);
    in32[1][0] = (uint32_t)(mask & In.y.re[0]);
    in32[2][0] = (uint32_t)(mask & In.z.re[0]);
    in32[3][0] = (uint32_t)(mask & In.t.re[0]);

    in32[0][1] = (uint32_t)(mask & ((In.x.re[0]>>28) + (In.x.re[1]<<27)));
    in32[1][1] = (uint32_t)(mask & ((In.y.re[0]>>28) + (In.y.re[1]<<27)));
    in32[2][1] = (uint32_t)(mask & ((In.z.re[0]>>28) + (In.z.re[1]<<27)));
    in32[3][1] = (uint32_t)(mask & ((In.t.re[0]>>28) + (In.t.re[1]<<27)));

    in32[0][2] = (uint32_t)(mask & (In.x.re[1]>>1));
    in32[1][2] = (uint32_t)(mask & (In.y.re[1]>>1));
    in32[2][2] = (uint32_t)(mask & (In.z.re[1]>>1));
    in32[3][2] = (uint32_t)(mask & (In.t.re[1]>>1));

    in32[0][3] = (uint32_t)(mask & ((In.x.re[1]>>29) + (In.x.re[2]<<26)));
    in32[1][3] = (uint32_t)(mask & ((In.y.re[1]>>29) + (In.y.re[2]<<26)));
    in32[2][3] = (uint32_t)(mask & ((In.z.re[1]>>29) + (In.z.re[2]<<26)));
    in32[3][3] = (uint32_t)(mask & ((In.t.re[1]>>29) + (In.t.re[2]<<26)));

    in32[0][4] = (uint32_t)(mask & (In.x.re[2]>>2));
    in32[1][4] = (uint32_t)(mask & (In.y.re[2]>>2));
    in32[2][4] = (uint32_t)(mask & (In.z.re[2]>>2));
    in32[3][4] = (uint32_t)(mask & (In.t.re[2]>>2));

    in32[0][5] = (uint32_t)(mask & ((In.x.re[2]>>30) + (In.x.re[3]<<25)));
    in32[1][5] = (uint32_t)(mask & ((In.y.re[2]>>30) + (In.y.re[3]<<25)));
    in32[2][5] = (uint32_t)(mask & ((In.z.re[2]>>30) + (In.z.re[3]<<25)));
    in32[3][5] = (uint32_t)(mask & ((In.t.re[2]>>30) + (In.t.re[3]<<25)));

    in32[0][6] = (uint32_t)(mask & (In.x.re[3]>>3));
    in32[1][6] = (uint32_t)(mask & (In.y.re[3]>>3));
    in32[2][6] = (uint32_t)(mask & (In.z.re[3]>>3));
    in32[3][6] = (uint32_t)(mask & (In.t.re[3]>>3));

    in32[0][7] = (uint32_t)(mask & ((In.x.re[3]>>31) + (In.x.re[4]<<24)));
    in32[1][7] = (uint32_t)(mask & ((In.y.re[3]>>31) + (In.y.re[4]<<24)));
    in32[2][7] = (uint32_t)(mask & ((In.z.re[3]>>31) + (In.z.re[4]<<24)));
    in32[3][7] = (uint32_t)(mask & ((In.t.re[3]>>31) + (In.t.re[4]<<24)));

    in32[0][8] = (uint32_t)(mask & (In.x.re[4]>>4));
    in32[1][8] = (uint32_t)(mask & (In.y.re[4]>>4));
    in32[2][8] = (uint32_t)(mask & (In.z.re[4]>>4));
    in32[3][8] = (uint32_t)(mask & (In.t.re[4]>>4));

    in32[0][9] = (uint32_t)(mask & ((In.x.re[4]>>32) + (In.x.re[5]<<23)));
    in32[1][9] = (uint32_t)(mask & ((In.y.re[4]>>32) + (In.y.re[5]<<23)));
    in32[2][9] = (uint32_t)(mask & ((In.z.re[4]>>32) + (In.z.re[5]<<23)));
    in32[3][9] = (uint32_t)(mask & ((In.t.re[4]>>32) + (In.t.re[5]<<23)));

    in32[0][10] = (uint32_t)(mask & (In.x.re[5]>>5));
    in32[1][10] = (uint32_t)(mask & (In.y.re[5]>>5));
    in32[2][10] = (uint32_t)(mask & (In.z.re[5]>>5));
    in32[3][10] = (uint32_t)(mask & (In.t.re[5]>>5));

    in32[0][11] = (uint32_t)(mask & ((In.x.re[5]>>33) + (In.x.re[6]<<22)));
    in32[1][11] = (uint32_t)(mask & ((In.y.re[5]>>33) + (In.y.re[6]<<22)));
    in32[2][11] = (uint32_t)(mask & ((In.z.re[5]>>33) + (In.z.re[6]<<22)));
    in32[3][11] = (uint32_t)(mask & ((In.t.re[5]>>33) + (In.t.re[6]<<22)));

    in32[0][12] = (uint32_t)(mask & (In.x.re[6]>>6));
    in32[1][12] = (uint32_t)(mask & (In.y.re[6]>>6));
    in32[2][12] = (uint32_t)(mask & (In.z.re[6]>>6));
    in32[3][12] = (uint32_t)(mask & (In.t.re[6]>>6));
    
    in32[0][13] = (uint32_t)(mask & (In.x.re[6]>>34));
    in32[1][13] = (uint32_t)(mask & (In.y.re[6]>>34));
    in32[2][13] = (uint32_t)(mask & (In.z.re[6]>>34));
    in32[3][13] = (uint32_t)(mask & (In.t.re[6]>>34));

    // im-cast
    in32[0][14] = (uint32_t)(mask & In.x.im[0]);
    in32[1][14] = (uint32_t)(mask & In.y.im[0]);
    in32[2][14] = (uint32_t)(mask & In.z.im[0]);
    in32[3][14] = (uint32_t)(mask & In.t.im[0]);

    in32[0][15] = (uint32_t)(mask & ((In.x.im[0]>>28) + (In.x.im[1]<<27)));
    in32[1][15] = (uint32_t)(mask & ((In.y.im[0]>>28) + (In.y.im[1]<<27)));
    in32[2][15] = (uint32_t)(mask & ((In.z.im[0]>>28) + (In.z.im[1]<<27)));
    in32[3][15] = (uint32_t)(mask & ((In.t.im[0]>>28) + (In.t.im[1]<<27)));

    in32[0][16] = (uint32_t)(mask & (In.x.im[1]>>1));
    in32[1][16] = (uint32_t)(mask & (In.y.im[1]>>1));
    in32[2][16] = (uint32_t)(mask & (In.z.im[1]>>1));
    in32[3][16] = (uint32_t)(mask & (In.t.im[1]>>1));

    in32[0][17] = (uint32_t)(mask & ((In.x.im[1]>>29) + (In.x.im[2]<<26)));
    in32[1][17] = (uint32_t)(mask & ((In.y.im[1]>>29) + (In.y.im[2]<<26)));
    in32[2][17] = (uint32_t)(mask & ((In.z.im[1]>>29) + (In.z.im[2]<<26)));
    in32[3][17] = (uint32_t)(mask & ((In.t.im[1]>>29) + (In.t.im[2]<<26)));

    in32[0][18] = (uint32_t)(mask & (In.x.im[2]>>2));
    in32[1][18] = (uint32_t)(mask & (In.y.im[2]>>2));
    in32[2][18] = (uint32_t)(mask & (In.z.im[2]>>2));
    in32[3][18] = (uint32_t)(mask & (In.t.im[2]>>2));

    in32[0][19] = (uint32_t)(mask & ((In.x.im[2]>>30) + (In.x.im[3]<<25)));
    in32[1][19] = (uint32_t)(mask & ((In.y.im[2]>>30) + (In.y.im[3]<<25)));
    in32[2][19] = (uint32_t)(mask & ((In.z.im[2]>>30) + (In.z.im[3]<<25)));
    in32[3][19] = (uint32_t)(mask & ((In.t.im[2]>>30) + (In.t.im[3]<<25)));

    in32[0][20] = (uint32_t)(mask & (In.x.im[3]>>3));
    in32[1][20] = (uint32_t)(mask & (In.y.im[3]>>3));
    in32[2][20] = (uint32_t)(mask & (In.z.im[3]>>3));
    in32[3][20] = (uint32_t)(mask & (In.t.im[3]>>3));

    in32[0][21] = (uint32_t)(mask & ((In.x.im[3]>>31) + (In.x.im[4]<<24)));
    in32[1][21] = (uint32_t)(mask & ((In.y.im[3]>>31) + (In.y.im[4]<<24)));
    in32[2][21] = (uint32_t)(mask & ((In.z.im[3]>>31) + (In.z.im[4]<<24)));
    in32[3][21] = (uint32_t)(mask & ((In.t.im[3]>>31) + (In.t.im[4]<<24)));

    in32[0][22] = (uint32_t)(mask & (In.x.im[4]>>4));
    in32[1][22] = (uint32_t)(mask & (In.y.im[4]>>4));
    in32[2][22] = (uint32_t)(mask & (In.z.im[4]>>4));
    in32[3][22] = (uint32_t)(mask & (In.t.im[4]>>4));

    in32[0][23] = (uint32_t)(mask & ((In.x.im[4]>>32) + (In.x.im[5]<<23)));
    in32[1][23] = (uint32_t)(mask & ((In.y.im[4]>>32) + (In.y.im[5]<<23)));
    in32[2][23] = (uint32_t)(mask & ((In.z.im[4]>>32) + (In.z.im[5]<<23)));
    in32[3][23] = (uint32_t)(mask & ((In.t.im[4]>>32) + (In.t.im[5]<<23)));

    in32[0][24] = (uint32_t)(mask & (In.x.im[5]>>5));
    in32[1][24] = (uint32_t)(mask & (In.y.im[5]>>5));
    in32[2][24] = (uint32_t)(mask & (In.z.im[5]>>5));
    in32[3][24] = (uint32_t)(mask & (In.t.im[5]>>5));

    in32[0][25] = (uint32_t)(mask & ((In.x.im[5]>>33) + (In.x.im[6]<<22)));
    in32[1][25] = (uint32_t)(mask & ((In.y.im[5]>>33) + (In.y.im[6]<<22)));
    in32[2][25] = (uint32_t)(mask & ((In.z.im[5]>>33) + (In.z.im[6]<<22)));
    in32[3][25] = (uint32_t)(mask & ((In.t.im[5]>>33) + (In.t.im[6]<<22)));

    in32[0][26] = (uint32_t)(mask & (In.x.im[6]>>6));
    in32[1][26] = (uint32_t)(mask & (In.y.im[6]>>6));
    in32[2][26] = (uint32_t)(mask & (In.z.im[6]>>6));
    in32[3][26] = (uint32_t)(mask & (In.t.im[6]>>6));
    
    in32[0][27] = (uint32_t)(mask & (In.x.im[6]>>34));
    in32[1][27] = (uint32_t)(mask & (In.y.im[6]>>34));
    in32[2][27] = (uint32_t)(mask & (In.z.im[6]>>34));
    in32[3][27] = (uint32_t)(mask & (In.t.im[6]>>34));

    // transpose
    for(int i = 0;i<FP2_LIMBS; i++){
        uint32x4_t tmp = {in32[0][i], in32[1][i], in32[2][i], in32[3][i]};
        Out[i] = tmp;
    }
}

void itranspose(theta_point_t *Out, uint32x4_t *In){
    uint32_t in32[4][FP2_LIMBS];

    // itranspose
    for(int i = 0;i<4;i++){
        for(int j = 0;j<FP2_LIMBS;j++){
            in32[i][j] = In[j][i];
        }
    }

    // re-cast
    uint64_t mask = ((uint64_t)1<<55)-1;
    Out->x.re[0] = mask & (((uint64_t)in32[0][0])  + ((uint64_t)in32[0][1] <<28));
    Out->x.im[0] = mask & (((uint64_t)in32[0][14]) + ((uint64_t)in32[0][15]<<28));
    Out->y.re[0] = mask & (((uint64_t)in32[1][0])  + ((uint64_t)in32[1][1] <<28));
    Out->y.im[0] = mask & (((uint64_t)in32[1][14]) + ((uint64_t)in32[1][15]<<28));
    Out->z.re[0] = mask & (((uint64_t)in32[2][0])  + ((uint64_t)in32[2][1] <<28));
    Out->z.im[0] = mask & (((uint64_t)in32[2][14]) + ((uint64_t)in32[2][15]<<28));
    Out->t.re[0] = mask & (((uint64_t)in32[3][0])  + ((uint64_t)in32[3][1] <<28));
    Out->t.im[0] = mask & (((uint64_t)in32[3][14]) + ((uint64_t)in32[3][15]<<28));

    Out->x.re[1] = mask & (((uint64_t)in32[0][1] >>27) + ((uint64_t)in32[0][2] <<1) + ((uint64_t)in32[0][3] <<29));
    Out->x.im[1] = mask & (((uint64_t)in32[0][15]>>27) + ((uint64_t)in32[0][16]<<1) + ((uint64_t)in32[0][17]<<29));
    Out->y.re[1] = mask & (((uint64_t)in32[1][1] >>27) + ((uint64_t)in32[1][2] <<1) + ((uint64_t)in32[1][3] <<29));
    Out->y.im[1] = mask & (((uint64_t)in32[1][15]>>27) + ((uint64_t)in32[1][16]<<1) + ((uint64_t)in32[1][17]<<29));
    Out->z.re[1] = mask & (((uint64_t)in32[2][1] >>27) + ((uint64_t)in32[2][2] <<1) + ((uint64_t)in32[2][3] <<29));
    Out->z.im[1] = mask & (((uint64_t)in32[2][15]>>27) + ((uint64_t)in32[2][16]<<1) + ((uint64_t)in32[2][17]<<29));
    Out->t.re[1] = mask & (((uint64_t)in32[3][1] >>27) + ((uint64_t)in32[3][2] <<1) + ((uint64_t)in32[3][3] <<29));
    Out->t.im[1] = mask & (((uint64_t)in32[3][15]>>27) + ((uint64_t)in32[3][16]<<1) + ((uint64_t)in32[3][17]<<29));

    Out->x.re[2] = mask & (((uint64_t)in32[0][3] >>26) + ((uint64_t)in32[0][4] <<2) + ((uint64_t)in32[0][5] <<30));
    Out->x.im[2] = mask & (((uint64_t)in32[0][17]>>26) + ((uint64_t)in32[0][18]<<2) + ((uint64_t)in32[0][19]<<30));
    Out->y.re[2] = mask & (((uint64_t)in32[1][3] >>26) + ((uint64_t)in32[1][4] <<2) + ((uint64_t)in32[1][5] <<30));
    Out->y.im[2] = mask & (((uint64_t)in32[1][17]>>26) + ((uint64_t)in32[1][18]<<2) + ((uint64_t)in32[1][19]<<30));
    Out->z.re[2] = mask & (((uint64_t)in32[2][3] >>26) + ((uint64_t)in32[2][4] <<2) + ((uint64_t)in32[2][5] <<30));
    Out->z.im[2] = mask & (((uint64_t)in32[2][17]>>26) + ((uint64_t)in32[2][18]<<2) + ((uint64_t)in32[2][19]<<30));
    Out->t.re[2] = mask & (((uint64_t)in32[3][3] >>26) + ((uint64_t)in32[3][4] <<2) + ((uint64_t)in32[3][5] <<30));
    Out->t.im[2] = mask & (((uint64_t)in32[3][17]>>26) + ((uint64_t)in32[3][18]<<2) + ((uint64_t)in32[3][19]<<30));

    Out->x.re[3] = mask & (((uint64_t)in32[0][5] >>25) + ((uint64_t)in32[0][6] <<3) + ((uint64_t)in32[0][7] <<31));
    Out->x.im[3] = mask & (((uint64_t)in32[0][19]>>25) + ((uint64_t)in32[0][20]<<3) + ((uint64_t)in32[0][21]<<31));
    Out->y.re[3] = mask & (((uint64_t)in32[1][5] >>25) + ((uint64_t)in32[1][6] <<3) + ((uint64_t)in32[1][7] <<31));
    Out->y.im[3] = mask & (((uint64_t)in32[1][19]>>25) + ((uint64_t)in32[1][20]<<3) + ((uint64_t)in32[1][21]<<31));
    Out->z.re[3] = mask & (((uint64_t)in32[2][5] >>25) + ((uint64_t)in32[2][6] <<3) + ((uint64_t)in32[2][7] <<31));
    Out->z.im[3] = mask & (((uint64_t)in32[2][19]>>25) + ((uint64_t)in32[2][20]<<3) + ((uint64_t)in32[2][21]<<31));
    Out->t.re[3] = mask & (((uint64_t)in32[3][5] >>25) + ((uint64_t)in32[3][6] <<3) + ((uint64_t)in32[3][7] <<31));
    Out->t.im[3] = mask & (((uint64_t)in32[3][19]>>25) + ((uint64_t)in32[3][20]<<3) + ((uint64_t)in32[3][21]<<31));

    Out->x.re[4] = mask & (((uint64_t)in32[0][7] >>24) + ((uint64_t)in32[0][8] <<4) + ((uint64_t)in32[0][9] <<32));
    Out->x.im[4] = mask & (((uint64_t)in32[0][21]>>24) + ((uint64_t)in32[0][22]<<4) + ((uint64_t)in32[0][23]<<32));
    Out->y.re[4] = mask & (((uint64_t)in32[1][7] >>24) + ((uint64_t)in32[1][8] <<4) + ((uint64_t)in32[1][9] <<32));
    Out->y.im[4] = mask & (((uint64_t)in32[1][21]>>24) + ((uint64_t)in32[1][22]<<4) + ((uint64_t)in32[1][23]<<32));
    Out->z.re[4] = mask & (((uint64_t)in32[2][7] >>24) + ((uint64_t)in32[2][8] <<4) + ((uint64_t)in32[2][9] <<32));
    Out->z.im[4] = mask & (((uint64_t)in32[2][21]>>24) + ((uint64_t)in32[2][22]<<4) + ((uint64_t)in32[2][23]<<32));
    Out->t.re[4] = mask & (((uint64_t)in32[3][7] >>24) + ((uint64_t)in32[3][8] <<4) + ((uint64_t)in32[3][9] <<32));
    Out->t.im[4] = mask & (((uint64_t)in32[3][21]>>24) + ((uint64_t)in32[3][22]<<4) + ((uint64_t)in32[3][23]<<32));

    Out->x.re[5] = mask & (((uint64_t)in32[0][9] >>23) + ((uint64_t)in32[0][10]<<5) + ((uint64_t)in32[0][11] <<33));
    Out->x.im[5] = mask & (((uint64_t)in32[0][23]>>23) + ((uint64_t)in32[0][24]<<5) + ((uint64_t)in32[0][25] <<33));
    Out->y.re[5] = mask & (((uint64_t)in32[1][9] >>23) + ((uint64_t)in32[1][10]<<5) + ((uint64_t)in32[1][11] <<33));
    Out->y.im[5] = mask & (((uint64_t)in32[1][23]>>23) + ((uint64_t)in32[1][24]<<5) + ((uint64_t)in32[1][25] <<33));
    Out->z.re[5] = mask & (((uint64_t)in32[2][9] >>23) + ((uint64_t)in32[2][10]<<5) + ((uint64_t)in32[2][11] <<33));
    Out->z.im[5] = mask & (((uint64_t)in32[2][23]>>23) + ((uint64_t)in32[2][24]<<5) + ((uint64_t)in32[2][25] <<33));
    Out->t.re[5] = mask & (((uint64_t)in32[3][9] >>23) + ((uint64_t)in32[3][10]<<5) + ((uint64_t)in32[3][11] <<33));
    Out->t.im[5] = mask & (((uint64_t)in32[3][23]>>23) + ((uint64_t)in32[3][24]<<5) + ((uint64_t)in32[3][25] <<33));

    Out->x.re[6] = mask & (((uint64_t)in32[0][11] >>22) + ((uint64_t)in32[0][12] <<6) + ((uint64_t)in32[0][13] <<34));
    Out->x.im[6] = mask & (((uint64_t)in32[0][25] >>22) + ((uint64_t)in32[0][26] <<6) + ((uint64_t)in32[0][27] <<34));
    Out->y.re[6] = mask & (((uint64_t)in32[1][11] >>22) + ((uint64_t)in32[1][12] <<6) + ((uint64_t)in32[1][13] <<34));
    Out->y.im[6] = mask & (((uint64_t)in32[1][25] >>22) + ((uint64_t)in32[1][26] <<6) + ((uint64_t)in32[1][27] <<34));
    Out->z.re[6] = mask & (((uint64_t)in32[2][11] >>22) + ((uint64_t)in32[2][12] <<6) + ((uint64_t)in32[2][13] <<34));
    Out->z.im[6] = mask & (((uint64_t)in32[2][25] >>22) + ((uint64_t)in32[2][26] <<6) + ((uint64_t)in32[2][27] <<34));
    Out->t.re[6] = mask & (((uint64_t)in32[3][11] >>22) + ((uint64_t)in32[3][12] <<6) + ((uint64_t)in32[3][13] <<34));
    Out->t.im[6] = mask & (((uint64_t)in32[3][25] >>22) + ((uint64_t)in32[3][26] <<6) + ((uint64_t)in32[3][27] <<34));
}