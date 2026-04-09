#include <stdint.h>
#include <stdio.h>
#include <arm_neon.h>
#include <hd.h>
#include <string.h>

void initial_montback_array(fp_t* mb, fp_t* imb){
    //2^513
    mb[0][0] = 303;
    mb[0][1] = 0;
    mb[0][2] = 0;
    mb[0][3] = 0;
    mb[0][4] = 0;
    mb[0][5] = 0;
    mb[0][6] = 0;
    mb[0][7] = 0;
    mb[0][8] = 193514046488576;
    //2^513
    imb[0][0] = 303;
    imb[0][1] = 0;
    imb[0][2] = 0;
    imb[0][3] = 0;
    imb[0][4] = 0;
    imb[0][5] = 0;
    imb[0][6] = 0;
    imb[0][7] = 0;
    imb[0][8] = 193514046488576;
}

void initial_q_value_1(uint32_t* q1){
    uint32_t tmp[FP_LIMBS] = {303, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 180224};
    memmove(q1, tmp, sizeof(uint32_t)*FP_LIMBS);
}

void ec_montback_array(fp_t* mb){
    //2^(513*4 - 513*3)
    mb[0][0] = 303;
    mb[0][1] = 0;
    mb[0][2] = 0;
    mb[0][3] = 0;
    mb[0][4] = 0;
    mb[0][5] = 0;
    mb[0][6] = 0;
    mb[0][7] = 0;
    mb[0][8] = 193514046488576;
}

void transpose(uint32x4_t *Out, theta_point_t In){
    uint64_t mask = ((uint64_t)1<<PER_LIMB)-1;
    uint32_t in32[4][FP2_LIMBS] = {0};

    // re-cast
    in32[0][0] = (uint32_t)(mask & In.x.re[0]);
    in32[1][0] = (uint32_t)(mask & In.y.re[0]);
    in32[2][0] = (uint32_t)(mask & In.z.re[0]);
    in32[3][0] = (uint32_t)(mask & In.t.re[0]);

    in32[0][1] = (uint32_t)(mask & ((In.x.re[0]>>27)));
    in32[1][1] = (uint32_t)(mask & ((In.y.re[0]>>27)));
    in32[2][1] = (uint32_t)(mask & ((In.z.re[0]>>27)));
    in32[3][1] = (uint32_t)(mask & ((In.t.re[0]>>27)));

    in32[0][2] = (uint32_t)(mask & ((In.x.re[0]>>54) + (In.x.re[1]<<3)));
    in32[1][2] = (uint32_t)(mask & ((In.y.re[0]>>54) + (In.y.re[1]<<3)));
    in32[2][2] = (uint32_t)(mask & ((In.z.re[0]>>54) + (In.z.re[1]<<3)));
    in32[3][2] = (uint32_t)(mask & ((In.t.re[0]>>54) + (In.t.re[1]<<3)));

    in32[0][3] = (uint32_t)(mask & (In.x.re[1]>>24));
    in32[1][3] = (uint32_t)(mask & (In.y.re[1]>>24));
    in32[2][3] = (uint32_t)(mask & (In.z.re[1]>>24));
    in32[3][3] = (uint32_t)(mask & (In.t.re[1]>>24));

    in32[0][4] = (uint32_t)(mask & ((In.x.re[1]>>51) + (In.x.re[2]<<6)));
    in32[1][4] = (uint32_t)(mask & ((In.y.re[1]>>51) + (In.y.re[2]<<6)));
    in32[2][4] = (uint32_t)(mask & ((In.z.re[1]>>51) + (In.z.re[2]<<6)));
    in32[3][4] = (uint32_t)(mask & ((In.t.re[1]>>51) + (In.t.re[2]<<6)));

    in32[0][5] = (uint32_t)(mask & (In.x.re[2]>>21));
    in32[1][5] = (uint32_t)(mask & (In.y.re[2]>>21));
    in32[2][5] = (uint32_t)(mask & (In.z.re[2]>>21));
    in32[3][5] = (uint32_t)(mask & (In.t.re[2]>>21));

    in32[0][6] = (uint32_t)(mask & ((In.x.re[2]>>48) + (In.x.re[3]<<9)));
    in32[1][6] = (uint32_t)(mask & ((In.y.re[2]>>48) + (In.y.re[3]<<9)));
    in32[2][6] = (uint32_t)(mask & ((In.z.re[2]>>48) + (In.z.re[3]<<9)));
    in32[3][6] = (uint32_t)(mask & ((In.t.re[2]>>48) + (In.t.re[3]<<9)));

    in32[0][7] = (uint32_t)(mask & (In.x.re[3]>>18));
    in32[1][7] = (uint32_t)(mask & (In.y.re[3]>>18));
    in32[2][7] = (uint32_t)(mask & (In.z.re[3]>>18));
    in32[3][7] = (uint32_t)(mask & (In.t.re[3]>>18));

    in32[0][8] = (uint32_t)(mask & ((In.x.re[3]>>45) + (In.x.re[4]<<12)));
    in32[1][8] = (uint32_t)(mask & ((In.y.re[3]>>45) + (In.y.re[4]<<12)));
    in32[2][8] = (uint32_t)(mask & ((In.z.re[3]>>45) + (In.z.re[4]<<12)));
    in32[3][8] = (uint32_t)(mask & ((In.t.re[3]>>45) + (In.t.re[4]<<12)));

    in32[0][9] = (uint32_t)(mask & (In.x.re[4]>>15));
    in32[1][9] = (uint32_t)(mask & (In.y.re[4]>>15));
    in32[2][9] = (uint32_t)(mask & (In.z.re[4]>>15));
    in32[3][9] = (uint32_t)(mask & (In.t.re[4]>>15));

    in32[0][10] = (uint32_t)(mask & ((In.x.re[4]>>42) + (In.x.re[5]<<15)));
    in32[1][10] = (uint32_t)(mask & ((In.y.re[4]>>42) + (In.y.re[5]<<15)));
    in32[2][10] = (uint32_t)(mask & ((In.z.re[4]>>42) + (In.z.re[5]<<15)));
    in32[3][10] = (uint32_t)(mask & ((In.t.re[4]>>42) + (In.t.re[5]<<15)));

    in32[0][11] = (uint32_t)(mask & (In.x.re[5]>>12));
    in32[1][11] = (uint32_t)(mask & (In.y.re[5]>>12));
    in32[2][11] = (uint32_t)(mask & (In.z.re[5]>>12));
    in32[3][11] = (uint32_t)(mask & (In.t.re[5]>>12));

    in32[0][12] = (uint32_t)(mask & + ((In.x.re[5]>>39) + (In.x.re[6]<<18)));
    in32[1][12] = (uint32_t)(mask & + ((In.y.re[5]>>39) + (In.y.re[6]<<18)));
    in32[2][12] = (uint32_t)(mask & + ((In.z.re[5]>>39) + (In.z.re[6]<<18)));
    in32[3][12] = (uint32_t)(mask & + ((In.t.re[5]>>39) + (In.t.re[6]<<18)));
    
    in32[0][13] = (uint32_t)(mask & (In.x.re[6]>>9));
    in32[1][13] = (uint32_t)(mask & (In.y.re[6]>>9));
    in32[2][13] = (uint32_t)(mask & (In.z.re[6]>>9));
    in32[3][13] = (uint32_t)(mask & (In.t.re[6]>>9));

    in32[0][14] = (uint32_t)(mask & ((In.x.re[6]>>36) + (In.x.re[7]<<21)));
    in32[1][14] = (uint32_t)(mask & ((In.y.re[6]>>36) + (In.y.re[7]<<21)));
    in32[2][14] = (uint32_t)(mask & ((In.z.re[6]>>36) + (In.z.re[7]<<21)));
    in32[3][14] = (uint32_t)(mask & ((In.t.re[6]>>36) + (In.t.re[7]<<21)));
    
    in32[0][15] = (uint32_t)(mask & (In.x.re[7]>>6));
    in32[1][15] = (uint32_t)(mask & (In.y.re[7]>>6));
    in32[2][15] = (uint32_t)(mask & (In.z.re[7]>>6));
    in32[3][15] = (uint32_t)(mask & (In.t.re[7]>>6));

    in32[0][16] = (uint32_t)(mask & ((In.x.re[7]>>33) + (In.x.re[8]<<24)));
    in32[1][16] = (uint32_t)(mask & ((In.y.re[7]>>33) + (In.y.re[8]<<24)));
    in32[2][16] = (uint32_t)(mask & ((In.z.re[7]>>33) + (In.z.re[8]<<24)));
    in32[3][16] = (uint32_t)(mask & ((In.t.re[7]>>33) + (In.t.re[8]<<24)));
    
    in32[0][17] = (uint32_t)(mask & (In.x.re[8]>>3));
    in32[1][17] = (uint32_t)(mask & (In.y.re[8]>>3));
    in32[2][17] = (uint32_t)(mask & (In.z.re[8]>>3));
    in32[3][17] = (uint32_t)(mask & (In.t.re[8]>>3));

    in32[0][18] = (uint32_t)(mask & (In.x.re[8]>>30));
    in32[1][18] = (uint32_t)(mask & (In.y.re[8]>>30));
    in32[2][18] = (uint32_t)(mask & (In.z.re[8]>>30));
    in32[3][18] = (uint32_t)(mask & (In.t.re[8]>>30));

    // im-cast
    in32[0][19] = (uint32_t)(mask & In.x.im[0]);
    in32[1][19] = (uint32_t)(mask & In.y.im[0]);
    in32[2][19] = (uint32_t)(mask & In.z.im[0]);
    in32[3][19] = (uint32_t)(mask & In.t.im[0]);

    in32[0][20] = (uint32_t)(mask & ((In.x.im[0]>>27)));
    in32[1][20] = (uint32_t)(mask & ((In.y.im[0]>>27)));
    in32[2][20] = (uint32_t)(mask & ((In.z.im[0]>>27)));
    in32[3][20] = (uint32_t)(mask & ((In.t.im[0]>>27)));

    in32[0][21] = (uint32_t)(mask & ((In.x.im[0]>>54) + (In.x.im[1]<<3)));
    in32[1][21] = (uint32_t)(mask & ((In.y.im[0]>>54) + (In.y.im[1]<<3)));
    in32[2][21] = (uint32_t)(mask & ((In.z.im[0]>>54) + (In.z.im[1]<<3)));
    in32[3][21] = (uint32_t)(mask & ((In.t.im[0]>>54) + (In.t.im[1]<<3)));

    in32[0][22] = (uint32_t)(mask & (In.x.im[1]>>24));
    in32[1][22] = (uint32_t)(mask & (In.y.im[1]>>24));
    in32[2][22] = (uint32_t)(mask & (In.z.im[1]>>24));
    in32[3][22] = (uint32_t)(mask & (In.t.im[1]>>24));

    in32[0][23] = (uint32_t)(mask & ((In.x.im[1]>>51) + (In.x.im[2]<<6)));
    in32[1][23] = (uint32_t)(mask & ((In.y.im[1]>>51) + (In.y.im[2]<<6)));
    in32[2][23] = (uint32_t)(mask & ((In.z.im[1]>>51) + (In.z.im[2]<<6)));
    in32[3][23] = (uint32_t)(mask & ((In.t.im[1]>>51) + (In.t.im[2]<<6)));

    in32[0][24] = (uint32_t)(mask & (In.x.im[2]>>21));
    in32[1][24] = (uint32_t)(mask & (In.y.im[2]>>21));
    in32[2][24] = (uint32_t)(mask & (In.z.im[2]>>21));
    in32[3][24] = (uint32_t)(mask & (In.t.im[2]>>21));

    in32[0][25] = (uint32_t)(mask & ((In.x.im[2]>>48) + (In.x.im[3]<<9)));
    in32[1][25] = (uint32_t)(mask & ((In.y.im[2]>>48) + (In.y.im[3]<<9)));
    in32[2][25] = (uint32_t)(mask & ((In.z.im[2]>>48) + (In.z.im[3]<<9)));
    in32[3][25] = (uint32_t)(mask & ((In.t.im[2]>>48) + (In.t.im[3]<<9)));

    in32[0][26] = (uint32_t)(mask & (In.x.im[3]>>18));
    in32[1][26] = (uint32_t)(mask & (In.y.im[3]>>18));
    in32[2][26] = (uint32_t)(mask & (In.z.im[3]>>18));
    in32[3][26] = (uint32_t)(mask & (In.t.im[3]>>18));

    in32[0][27] = (uint32_t)(mask & ((In.x.im[3]>>45) + (In.x.im[4]<<12)));
    in32[1][27] = (uint32_t)(mask & ((In.y.im[3]>>45) + (In.y.im[4]<<12)));
    in32[2][27] = (uint32_t)(mask & ((In.z.im[3]>>45) + (In.z.im[4]<<12)));
    in32[3][27] = (uint32_t)(mask & ((In.t.im[3]>>45) + (In.t.im[4]<<12)));

    in32[0][28] = (uint32_t)(mask & (In.x.im[4]>>15));
    in32[1][28] = (uint32_t)(mask & (In.y.im[4]>>15));
    in32[2][28] = (uint32_t)(mask & (In.z.im[4]>>15));
    in32[3][28] = (uint32_t)(mask & (In.t.im[4]>>15));

    in32[0][29] = (uint32_t)(mask & ((In.x.im[4]>>42) + (In.x.im[5]<<15)));
    in32[1][29] = (uint32_t)(mask & ((In.y.im[4]>>42) + (In.y.im[5]<<15)));
    in32[2][29] = (uint32_t)(mask & ((In.z.im[4]>>42) + (In.z.im[5]<<15)));
    in32[3][29] = (uint32_t)(mask & ((In.t.im[4]>>42) + (In.t.im[5]<<15)));

    in32[0][30] = (uint32_t)(mask & (In.x.im[5]>>12));
    in32[1][30] = (uint32_t)(mask & (In.y.im[5]>>12));
    in32[2][30] = (uint32_t)(mask & (In.z.im[5]>>12));
    in32[3][30] = (uint32_t)(mask & (In.t.im[5]>>12));

    in32[0][31] = (uint32_t)(mask & + ((In.x.im[5]>>39) + (In.x.im[6]<<18)));
    in32[1][31] = (uint32_t)(mask & + ((In.y.im[5]>>39) + (In.y.im[6]<<18)));
    in32[2][31] = (uint32_t)(mask & + ((In.z.im[5]>>39) + (In.z.im[6]<<18)));
    in32[3][31] = (uint32_t)(mask & + ((In.t.im[5]>>39) + (In.t.im[6]<<18)));
    
    in32[0][32] = (uint32_t)(mask & (In.x.im[6]>>9));
    in32[1][32] = (uint32_t)(mask & (In.y.im[6]>>9));
    in32[2][32] = (uint32_t)(mask & (In.z.im[6]>>9));
    in32[3][32] = (uint32_t)(mask & (In.t.im[6]>>9));

    in32[0][33] = (uint32_t)(mask & ((In.x.im[6]>>36) + (In.x.im[7]<<21)));
    in32[1][33] = (uint32_t)(mask & ((In.y.im[6]>>36) + (In.y.im[7]<<21)));
    in32[2][33] = (uint32_t)(mask & ((In.z.im[6]>>36) + (In.z.im[7]<<21)));
    in32[3][33] = (uint32_t)(mask & ((In.t.im[6]>>36) + (In.t.im[7]<<21)));
    
    in32[0][34] = (uint32_t)(mask & (In.x.im[7]>>6));
    in32[1][34] = (uint32_t)(mask & (In.y.im[7]>>6));
    in32[2][34] = (uint32_t)(mask & (In.z.im[7]>>6));
    in32[3][34] = (uint32_t)(mask & (In.t.im[7]>>6));

    in32[0][35] = (uint32_t)(mask & ((In.x.im[7]>>33) + (In.x.im[8]<<24)));
    in32[1][35] = (uint32_t)(mask & ((In.y.im[7]>>33) + (In.y.im[8]<<24)));
    in32[2][35] = (uint32_t)(mask & ((In.z.im[7]>>33) + (In.z.im[8]<<24)));
    in32[3][35] = (uint32_t)(mask & ((In.t.im[7]>>33) + (In.t.im[8]<<24)));
    
    in32[0][36] = (uint32_t)(mask & (In.x.im[8]>>3));
    in32[1][36] = (uint32_t)(mask & (In.y.im[8]>>3));
    in32[2][36] = (uint32_t)(mask & (In.z.im[8]>>3));
    in32[3][36] = (uint32_t)(mask & (In.t.im[8]>>3));

    in32[0][37] = (uint32_t)(mask & (In.x.im[8]>>30));
    in32[1][37] = (uint32_t)(mask & (In.y.im[8]>>30));
    in32[2][37] = (uint32_t)(mask & (In.z.im[8]>>30));
    in32[3][37] = (uint32_t)(mask & (In.t.im[8]>>30));

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

    uint64_t mask = ((uint64_t)1<<57)-1;
    Out->x.re[0] = mask & (((uint64_t)in32[0][0])  + ((uint64_t)in32[0][1] <<27) + ((uint64_t)in32[0][2] <<54));
    Out->x.im[0] = mask & (((uint64_t)in32[0][19]) + ((uint64_t)in32[0][20]<<27) + ((uint64_t)in32[0][21]<<54));
    Out->y.re[0] = mask & (((uint64_t)in32[1][0])  + ((uint64_t)in32[1][1] <<27) + ((uint64_t)in32[1][2] <<54));
    Out->y.im[0] = mask & (((uint64_t)in32[1][19]) + ((uint64_t)in32[1][20]<<27) + ((uint64_t)in32[1][21]<<54));
    Out->z.re[0] = mask & (((uint64_t)in32[2][0])  + ((uint64_t)in32[2][1] <<27) + ((uint64_t)in32[2][2] <<54));
    Out->z.im[0] = mask & (((uint64_t)in32[2][19]) + ((uint64_t)in32[2][20]<<27) + ((uint64_t)in32[2][21]<<54));
    Out->t.re[0] = mask & (((uint64_t)in32[3][0])  + ((uint64_t)in32[3][1] <<27) + ((uint64_t)in32[3][2] <<54));
    Out->t.im[0] = mask & (((uint64_t)in32[3][19]) + ((uint64_t)in32[3][20]<<27) + ((uint64_t)in32[3][21]<<54));

    Out->x.re[1] = mask & (((uint64_t)in32[0][2] >>3) + ((uint64_t)in32[0][3] <<24) + ((uint64_t)in32[0][4] <<51));
    Out->x.im[1] = mask & (((uint64_t)in32[0][21]>>3) + ((uint64_t)in32[0][22]<<24) + ((uint64_t)in32[0][23]<<51));
    Out->y.re[1] = mask & (((uint64_t)in32[1][2] >>3) + ((uint64_t)in32[1][3] <<24) + ((uint64_t)in32[1][4] <<51));
    Out->y.im[1] = mask & (((uint64_t)in32[1][21]>>3) + ((uint64_t)in32[1][22]<<24) + ((uint64_t)in32[1][23]<<51));
    Out->z.re[1] = mask & (((uint64_t)in32[2][2] >>3) + ((uint64_t)in32[2][3] <<24) + ((uint64_t)in32[2][4] <<51));
    Out->z.im[1] = mask & (((uint64_t)in32[2][21]>>3) + ((uint64_t)in32[2][22]<<24) + ((uint64_t)in32[2][23]<<51));
    Out->t.re[1] = mask & (((uint64_t)in32[3][2] >>3) + ((uint64_t)in32[3][3] <<24) + ((uint64_t)in32[3][4] <<51));
    Out->t.im[1] = mask & (((uint64_t)in32[3][21]>>3) + ((uint64_t)in32[3][22]<<24) + ((uint64_t)in32[3][23]<<51));

    Out->x.re[2] = mask & (((uint64_t)in32[0][4] >>6) + ((uint64_t)in32[0][5] <<21) + ((uint64_t)in32[0][6] <<48));
    Out->x.im[2] = mask & (((uint64_t)in32[0][23]>>6) + ((uint64_t)in32[0][24]<<21) + ((uint64_t)in32[0][25]<<48));
    Out->y.re[2] = mask & (((uint64_t)in32[1][4] >>6) + ((uint64_t)in32[1][5] <<21) + ((uint64_t)in32[1][6] <<48));
    Out->y.im[2] = mask & (((uint64_t)in32[1][23]>>6) + ((uint64_t)in32[1][24]<<21) + ((uint64_t)in32[1][25]<<48));
    Out->z.re[2] = mask & (((uint64_t)in32[2][4] >>6) + ((uint64_t)in32[2][5] <<21) + ((uint64_t)in32[2][6] <<48));
    Out->z.im[2] = mask & (((uint64_t)in32[2][23]>>6) + ((uint64_t)in32[2][24]<<21) + ((uint64_t)in32[2][25]<<48));
    Out->t.re[2] = mask & (((uint64_t)in32[3][4] >>6) + ((uint64_t)in32[3][5] <<21) + ((uint64_t)in32[3][6] <<48));
    Out->t.im[2] = mask & (((uint64_t)in32[3][23]>>6) + ((uint64_t)in32[3][24]<<21) + ((uint64_t)in32[3][25]<<48));

    Out->x.re[3] = mask & (((uint64_t)in32[0][6] >>9) + ((uint64_t)in32[0][7] <<18) + ((uint64_t)in32[0][8] <<45));
    Out->x.im[3] = mask & (((uint64_t)in32[0][25]>>9) + ((uint64_t)in32[0][26]<<18) + ((uint64_t)in32[0][27]<<45));
    Out->y.re[3] = mask & (((uint64_t)in32[1][6] >>9) + ((uint64_t)in32[1][7] <<18) + ((uint64_t)in32[1][8] <<45));
    Out->y.im[3] = mask & (((uint64_t)in32[1][25]>>9) + ((uint64_t)in32[1][26]<<18) + ((uint64_t)in32[1][27]<<45));
    Out->z.re[3] = mask & (((uint64_t)in32[2][6] >>9) + ((uint64_t)in32[2][7] <<18) + ((uint64_t)in32[2][8] <<45));
    Out->z.im[3] = mask & (((uint64_t)in32[2][25]>>9) + ((uint64_t)in32[2][26]<<18) + ((uint64_t)in32[2][27]<<45));
    Out->t.re[3] = mask & (((uint64_t)in32[3][6] >>9) + ((uint64_t)in32[3][7] <<18) + ((uint64_t)in32[3][8] <<45));
    Out->t.im[3] = mask & (((uint64_t)in32[3][25]>>9) + ((uint64_t)in32[3][26]<<18) + ((uint64_t)in32[3][27]<<45));

    Out->x.re[4] = mask & (((uint64_t)in32[0][8] >>12) + ((uint64_t)in32[0][9] <<15) + ((uint64_t)in32[0][10]<<42));
    Out->x.im[4] = mask & (((uint64_t)in32[0][27]>>12) + ((uint64_t)in32[0][28]<<15) + ((uint64_t)in32[0][29]<<42));
    Out->y.re[4] = mask & (((uint64_t)in32[1][8] >>12) + ((uint64_t)in32[1][9] <<15) + ((uint64_t)in32[1][10]<<42));
    Out->y.im[4] = mask & (((uint64_t)in32[1][27]>>12) + ((uint64_t)in32[1][28]<<15) + ((uint64_t)in32[1][29]<<42));
    Out->z.re[4] = mask & (((uint64_t)in32[2][8] >>12) + ((uint64_t)in32[2][9] <<15) + ((uint64_t)in32[2][10]<<42));
    Out->z.im[4] = mask & (((uint64_t)in32[2][27]>>12) + ((uint64_t)in32[2][28]<<15) + ((uint64_t)in32[2][29]<<42));
    Out->t.re[4] = mask & (((uint64_t)in32[3][8] >>12) + ((uint64_t)in32[3][9] <<15) + ((uint64_t)in32[3][10]<<42));
    Out->t.im[4] = mask & (((uint64_t)in32[3][27]>>12) + ((uint64_t)in32[3][28]<<15) + ((uint64_t)in32[3][29]<<42));

    Out->x.re[5] = mask & (((uint64_t)in32[0][10]>>15) + ((uint64_t)in32[0][11]<<12) + ((uint64_t)in32[0][12]<<39));
    Out->x.im[5] = mask & (((uint64_t)in32[0][29]>>15) + ((uint64_t)in32[0][30]<<12) + ((uint64_t)in32[0][31]<<39));
    Out->y.re[5] = mask & (((uint64_t)in32[1][10]>>15) + ((uint64_t)in32[1][11]<<12) + ((uint64_t)in32[1][12]<<39));
    Out->y.im[5] = mask & (((uint64_t)in32[1][29]>>15) + ((uint64_t)in32[1][30]<<12) + ((uint64_t)in32[1][31]<<39));
    Out->z.re[5] = mask & (((uint64_t)in32[2][10]>>15) + ((uint64_t)in32[2][11]<<12) + ((uint64_t)in32[2][12]<<39));
    Out->z.im[5] = mask & (((uint64_t)in32[2][29]>>15) + ((uint64_t)in32[2][30]<<12) + ((uint64_t)in32[2][31]<<39));
    Out->t.re[5] = mask & (((uint64_t)in32[3][10]>>15) + ((uint64_t)in32[3][11]<<12) + ((uint64_t)in32[3][12]<<39));
    Out->t.im[5] = mask & (((uint64_t)in32[3][29]>>15) + ((uint64_t)in32[3][30]<<12) + ((uint64_t)in32[3][31]<<39));

    Out->x.re[6] = mask & (((uint64_t)in32[0][12]>>18) + ((uint64_t)in32[0][13]<<9) + ((uint64_t)in32[0][14]<<36));
    Out->x.im[6] = mask & (((uint64_t)in32[0][31]>>18) + ((uint64_t)in32[0][32]<<9) + ((uint64_t)in32[0][33]<<36));
    Out->y.re[6] = mask & (((uint64_t)in32[1][12]>>18) + ((uint64_t)in32[1][13]<<9) + ((uint64_t)in32[1][14]<<36));
    Out->y.im[6] = mask & (((uint64_t)in32[1][31]>>18) + ((uint64_t)in32[1][32]<<9) + ((uint64_t)in32[1][33]<<36));
    Out->z.re[6] = mask & (((uint64_t)in32[2][12]>>18) + ((uint64_t)in32[2][13]<<9) + ((uint64_t)in32[2][14]<<36));
    Out->z.im[6] = mask & (((uint64_t)in32[2][31]>>18) + ((uint64_t)in32[2][32]<<9) + ((uint64_t)in32[2][33]<<36));
    Out->t.re[6] = mask & (((uint64_t)in32[3][12]>>18) + ((uint64_t)in32[3][13]<<9) + ((uint64_t)in32[3][14]<<36));
    Out->t.im[6] = mask & (((uint64_t)in32[3][31]>>18) + ((uint64_t)in32[3][32]<<9) + ((uint64_t)in32[3][33]<<36));

    Out->x.re[7] = mask & (((uint64_t)in32[0][14]>>21) + ((uint64_t)in32[0][15]<<6) + ((uint64_t)in32[0][16]<<33));
    Out->x.im[7] = mask & (((uint64_t)in32[0][33]>>21) + ((uint64_t)in32[0][34]<<6) + ((uint64_t)in32[0][35]<<33));
    Out->y.re[7] = mask & (((uint64_t)in32[1][14]>>21) + ((uint64_t)in32[1][15]<<6) + ((uint64_t)in32[1][16]<<33));
    Out->y.im[7] = mask & (((uint64_t)in32[1][33]>>21) + ((uint64_t)in32[1][34]<<6) + ((uint64_t)in32[1][35]<<33));
    Out->z.re[7] = mask & (((uint64_t)in32[2][14]>>21) + ((uint64_t)in32[2][15]<<6) + ((uint64_t)in32[2][16]<<33));
    Out->z.im[7] = mask & (((uint64_t)in32[2][33]>>21) + ((uint64_t)in32[2][34]<<6) + ((uint64_t)in32[2][35]<<33));
    Out->t.re[7] = mask & (((uint64_t)in32[3][14]>>21) + ((uint64_t)in32[3][15]<<6) + ((uint64_t)in32[3][16]<<33));
    Out->t.im[7] = mask & (((uint64_t)in32[3][33]>>21) + ((uint64_t)in32[3][34]<<6) + ((uint64_t)in32[3][35]<<33));

    Out->x.re[8] = mask & (((uint64_t)in32[0][16]>>24) + ((uint64_t)in32[0][17]<<3) + ((uint64_t)in32[0][18]<<30));
    Out->x.im[8] = mask & (((uint64_t)in32[0][35]>>24) + ((uint64_t)in32[0][36]<<3) + ((uint64_t)in32[0][37]<<30));
    Out->y.re[8] = mask & (((uint64_t)in32[1][16]>>24) + ((uint64_t)in32[1][17]<<3) + ((uint64_t)in32[1][18]<<30));
    Out->y.im[8] = mask & (((uint64_t)in32[1][35]>>24) + ((uint64_t)in32[1][36]<<3) + ((uint64_t)in32[1][37]<<30));
    Out->z.re[8] = mask & (((uint64_t)in32[2][16]>>24) + ((uint64_t)in32[2][17]<<3) + ((uint64_t)in32[2][18]<<30));
    Out->z.im[8] = mask & (((uint64_t)in32[2][35]>>24) + ((uint64_t)in32[2][36]<<3) + ((uint64_t)in32[2][37]<<30));
    Out->t.re[8] = mask & (((uint64_t)in32[3][16]>>24) + ((uint64_t)in32[3][17]<<3) + ((uint64_t)in32[3][18]<<30));
    Out->t.im[8] = mask & (((uint64_t)in32[3][35]>>24) + ((uint64_t)in32[3][36]<<3) + ((uint64_t)in32[3][37]<<30));
}

