

#include "Types.h"
#include "BRIEF.h"
#include <cmath>
#include <iostream>
#include "immintrin.h"


extern "C" void kernel1(int *p_x_a, int *p_y_a, int *p_x_b, int *p_y_b,
                               int *i_x_a, int *i_y_a, int *i_x_b, int *i_y_b,
                               float cos, float sin);

void
BRIEF::rbrief(unsigned char *image_src, const int height_image, const int width_image, const int n_channels,
                const int stride_image, const int *x, const int *y, const float *angle, const int n_features, int64_t *bd,
                const int n_rows_bd) {
    for (int j = 0; j < n_features; j++) {
        if ((x[j] > diag_length_pattern) && x[j] < (width_image - diag_length_pattern)
            && (y[j] > diag_length_pattern) && y[j] < (height_image - diag_length_pattern)) {
            float cos_angle = std::cos(angle[j]);
            float sin_angle = std::sin(angle[j]);
            float tri0[4]  __attribute__((aligned(32))) = {cos_angle, sin_angle, cos_angle, sin_angle};
            float tri1[4]  __attribute__((aligned(32))) = {sin_angle,cos_angle, sin_angle, cos_angle};
            unsigned char *image_center = image_src + y[j] * stride_image + x[j] * n_channels;
            // N_DIM_BINARYDESCRIPTOR / SIZE_BITS_HAMING = 4
            unsigned char a[256] __attribute__((aligned(32)));
            unsigned char b[256] __attribute__((aligned(32)));
            int32_t i_x_a[256] __attribute__((aligned(32)));
            int32_t i_y_a[256] __attribute__((aligned(32)));
            int32_t i_x_b[256] __attribute__((aligned(32)));
            int32_t i_y_b[256] __attribute__((aligned(32)));
            //kernel1(gaussian_bit_pattern_31_x_a,gaussian_bit_pattern_31_y_a,gaussian_bit_pattern_31_x_b,gaussian_bit_pattern_31_y_b,
            //    i_x_a,i_y_a, i_x_b,i_y_b, cos_angle, sin_angle);
            for (int i = 0 ; i< 256; i++)  {
                    i_x_a[i] = (int)(gaussian_bit_pattern_31_x_a[i] * cos_angle - gaussian_bit_pattern_31_y_a[i] * sin_angle);
                    i_y_a[i] = (int)(gaussian_bit_pattern_31_x_a[i] * sin_angle + gaussian_bit_pattern_31_y_a[i] * cos_angle);
                    i_x_b[i] = (int)(gaussian_bit_pattern_31_x_b[i] * cos_angle - gaussian_bit_pattern_31_y_b[i] * sin_angle);
                    i_y_b[i] = (int)(gaussian_bit_pattern_31_x_b[i] * sin_angle + gaussian_bit_pattern_31_y_b[i] * cos_angle);
            }
            a[0]= *(image_center + i_y_a[0]*stride_image + i_x_a[0]*n_channels); b[0]= *(image_center + i_y_b[0]*stride_image + i_x_b[0]*n_channels);
            a[1]= *(image_center + i_y_a[1]*stride_image + i_x_a[1]*n_channels); b[1]= *(image_center + i_y_b[1]*stride_image + i_x_b[1]*n_channels);
            a[2]= *(image_center + i_y_a[2]*stride_image + i_x_a[2]*n_channels); b[2]= *(image_center + i_y_b[2]*stride_image + i_x_b[2]*n_channels);
            a[3]= *(image_center + i_y_a[3]*stride_image + i_x_a[3]*n_channels); b[3]= *(image_center + i_y_b[3]*stride_image + i_x_b[3]*n_channels);
            a[4]= *(image_center + i_y_a[4]*stride_image + i_x_a[4]*n_channels); b[4]= *(image_center + i_y_b[4]*stride_image + i_x_b[4]*n_channels);
            a[5]= *(image_center + i_y_a[5]*stride_image + i_x_a[5]*n_channels); b[5]= *(image_center + i_y_b[5]*stride_image + i_x_b[5]*n_channels);
            a[6]= *(image_center + i_y_a[6]*stride_image + i_x_a[6]*n_channels); b[6]= *(image_center + i_y_b[6]*stride_image + i_x_b[6]*n_channels);
            a[7]= *(image_center + i_y_a[7]*stride_image + i_x_a[7]*n_channels); b[7]= *(image_center + i_y_b[7]*stride_image + i_x_b[7]*n_channels);
            a[8]= *(image_center + i_y_a[8]*stride_image + i_x_a[8]*n_channels); b[8]= *(image_center + i_y_b[8]*stride_image + i_x_b[8]*n_channels);
            a[9]= *(image_center + i_y_a[9]*stride_image + i_x_a[9]*n_channels); b[9]= *(image_center + i_y_b[9]*stride_image + i_x_b[9]*n_channels);
            a[10]= *(image_center + i_y_a[10]*stride_image + i_x_a[10]*n_channels); b[10]= *(image_center + i_y_b[10]*stride_image + i_x_b[10]*n_channels);
            a[11]= *(image_center + i_y_a[11]*stride_image + i_x_a[11]*n_channels); b[11]= *(image_center + i_y_b[11]*stride_image + i_x_b[11]*n_channels);
            a[12]= *(image_center + i_y_a[12]*stride_image + i_x_a[12]*n_channels); b[12]= *(image_center + i_y_b[12]*stride_image + i_x_b[12]*n_channels);
            a[13]= *(image_center + i_y_a[13]*stride_image + i_x_a[13]*n_channels); b[13]= *(image_center + i_y_b[13]*stride_image + i_x_b[13]*n_channels);
            a[14]= *(image_center + i_y_a[14]*stride_image + i_x_a[14]*n_channels); b[14]= *(image_center + i_y_b[14]*stride_image + i_x_b[14]*n_channels);
            a[15]= *(image_center + i_y_a[15]*stride_image + i_x_a[15]*n_channels); b[15]= *(image_center + i_y_b[15]*stride_image + i_x_b[15]*n_channels);
            a[16]= *(image_center + i_y_a[16]*stride_image + i_x_a[16]*n_channels); b[16]= *(image_center + i_y_b[16]*stride_image + i_x_b[16]*n_channels);
            a[17]= *(image_center + i_y_a[17]*stride_image + i_x_a[17]*n_channels); b[17]= *(image_center + i_y_b[17]*stride_image + i_x_b[17]*n_channels);
            a[18]= *(image_center + i_y_a[18]*stride_image + i_x_a[18]*n_channels); b[18]= *(image_center + i_y_b[18]*stride_image + i_x_b[18]*n_channels);
            a[19]= *(image_center + i_y_a[19]*stride_image + i_x_a[19]*n_channels); b[19]= *(image_center + i_y_b[19]*stride_image + i_x_b[19]*n_channels);
            a[20]= *(image_center + i_y_a[20]*stride_image + i_x_a[20]*n_channels); b[20]= *(image_center + i_y_b[20]*stride_image + i_x_b[20]*n_channels);
            a[21]= *(image_center + i_y_a[21]*stride_image + i_x_a[21]*n_channels); b[21]= *(image_center + i_y_b[21]*stride_image + i_x_b[21]*n_channels);
            a[22]= *(image_center + i_y_a[22]*stride_image + i_x_a[22]*n_channels); b[22]= *(image_center + i_y_b[22]*stride_image + i_x_b[22]*n_channels);
            a[23]= *(image_center + i_y_a[23]*stride_image + i_x_a[23]*n_channels); b[23]= *(image_center + i_y_b[23]*stride_image + i_x_b[23]*n_channels);
            a[24]= *(image_center + i_y_a[24]*stride_image + i_x_a[24]*n_channels); b[24]= *(image_center + i_y_b[24]*stride_image + i_x_b[24]*n_channels);
            a[25]= *(image_center + i_y_a[25]*stride_image + i_x_a[25]*n_channels); b[25]= *(image_center + i_y_b[25]*stride_image + i_x_b[25]*n_channels);
            a[26]= *(image_center + i_y_a[26]*stride_image + i_x_a[26]*n_channels); b[26]= *(image_center + i_y_b[26]*stride_image + i_x_b[26]*n_channels);
            a[27]= *(image_center + i_y_a[27]*stride_image + i_x_a[27]*n_channels); b[27]= *(image_center + i_y_b[27]*stride_image + i_x_b[27]*n_channels);
            a[28]= *(image_center + i_y_a[28]*stride_image + i_x_a[28]*n_channels); b[28]= *(image_center + i_y_b[28]*stride_image + i_x_b[28]*n_channels);
            a[29]= *(image_center + i_y_a[29]*stride_image + i_x_a[29]*n_channels); b[29]= *(image_center + i_y_b[29]*stride_image + i_x_b[29]*n_channels);
            a[30]= *(image_center + i_y_a[30]*stride_image + i_x_a[30]*n_channels); b[30]= *(image_center + i_y_b[30]*stride_image + i_x_b[30]*n_channels);
            a[31]= *(image_center + i_y_a[31]*stride_image + i_x_a[31]*n_channels); b[31]= *(image_center + i_y_b[31]*stride_image + i_x_b[31]*n_channels);
            a[32]= *(image_center + i_y_a[32]*stride_image + i_x_a[32]*n_channels); b[32]= *(image_center + i_y_b[32]*stride_image + i_x_b[32]*n_channels);
            a[33]= *(image_center + i_y_a[33]*stride_image + i_x_a[33]*n_channels); b[33]= *(image_center + i_y_b[33]*stride_image + i_x_b[33]*n_channels);
            a[34]= *(image_center + i_y_a[34]*stride_image + i_x_a[34]*n_channels); b[34]= *(image_center + i_y_b[34]*stride_image + i_x_b[34]*n_channels);
            a[35]= *(image_center + i_y_a[35]*stride_image + i_x_a[35]*n_channels); b[35]= *(image_center + i_y_b[35]*stride_image + i_x_b[35]*n_channels);
            a[36]= *(image_center + i_y_a[36]*stride_image + i_x_a[36]*n_channels); b[36]= *(image_center + i_y_b[36]*stride_image + i_x_b[36]*n_channels);
            a[37]= *(image_center + i_y_a[37]*stride_image + i_x_a[37]*n_channels); b[37]= *(image_center + i_y_b[37]*stride_image + i_x_b[37]*n_channels);
            a[38]= *(image_center + i_y_a[38]*stride_image + i_x_a[38]*n_channels); b[38]= *(image_center + i_y_b[38]*stride_image + i_x_b[38]*n_channels);
            a[39]= *(image_center + i_y_a[39]*stride_image + i_x_a[39]*n_channels); b[39]= *(image_center + i_y_b[39]*stride_image + i_x_b[39]*n_channels);
            a[40]= *(image_center + i_y_a[40]*stride_image + i_x_a[40]*n_channels); b[40]= *(image_center + i_y_b[40]*stride_image + i_x_b[40]*n_channels);
            a[41]= *(image_center + i_y_a[41]*stride_image + i_x_a[41]*n_channels); b[41]= *(image_center + i_y_b[41]*stride_image + i_x_b[41]*n_channels);
            a[42]= *(image_center + i_y_a[42]*stride_image + i_x_a[42]*n_channels); b[42]= *(image_center + i_y_b[42]*stride_image + i_x_b[42]*n_channels);
            a[43]= *(image_center + i_y_a[43]*stride_image + i_x_a[43]*n_channels); b[43]= *(image_center + i_y_b[43]*stride_image + i_x_b[43]*n_channels);
            a[44]= *(image_center + i_y_a[44]*stride_image + i_x_a[44]*n_channels); b[44]= *(image_center + i_y_b[44]*stride_image + i_x_b[44]*n_channels);
            a[45]= *(image_center + i_y_a[45]*stride_image + i_x_a[45]*n_channels); b[45]= *(image_center + i_y_b[45]*stride_image + i_x_b[45]*n_channels);
            a[46]= *(image_center + i_y_a[46]*stride_image + i_x_a[46]*n_channels); b[46]= *(image_center + i_y_b[46]*stride_image + i_x_b[46]*n_channels);
            a[47]= *(image_center + i_y_a[47]*stride_image + i_x_a[47]*n_channels); b[47]= *(image_center + i_y_b[47]*stride_image + i_x_b[47]*n_channels);
            a[48]= *(image_center + i_y_a[48]*stride_image + i_x_a[48]*n_channels); b[48]= *(image_center + i_y_b[48]*stride_image + i_x_b[48]*n_channels);
            a[49]= *(image_center + i_y_a[49]*stride_image + i_x_a[49]*n_channels); b[49]= *(image_center + i_y_b[49]*stride_image + i_x_b[49]*n_channels);
            a[50]= *(image_center + i_y_a[50]*stride_image + i_x_a[50]*n_channels); b[50]= *(image_center + i_y_b[50]*stride_image + i_x_b[50]*n_channels);
            a[51]= *(image_center + i_y_a[51]*stride_image + i_x_a[51]*n_channels); b[51]= *(image_center + i_y_b[51]*stride_image + i_x_b[51]*n_channels);
            a[52]= *(image_center + i_y_a[52]*stride_image + i_x_a[52]*n_channels); b[52]= *(image_center + i_y_b[52]*stride_image + i_x_b[52]*n_channels);
            a[53]= *(image_center + i_y_a[53]*stride_image + i_x_a[53]*n_channels); b[53]= *(image_center + i_y_b[53]*stride_image + i_x_b[53]*n_channels);
            a[54]= *(image_center + i_y_a[54]*stride_image + i_x_a[54]*n_channels); b[54]= *(image_center + i_y_b[54]*stride_image + i_x_b[54]*n_channels);
            a[55]= *(image_center + i_y_a[55]*stride_image + i_x_a[55]*n_channels); b[55]= *(image_center + i_y_b[55]*stride_image + i_x_b[55]*n_channels);
            a[56]= *(image_center + i_y_a[56]*stride_image + i_x_a[56]*n_channels); b[56]= *(image_center + i_y_b[56]*stride_image + i_x_b[56]*n_channels);
            a[57]= *(image_center + i_y_a[57]*stride_image + i_x_a[57]*n_channels); b[57]= *(image_center + i_y_b[57]*stride_image + i_x_b[57]*n_channels);
            a[58]= *(image_center + i_y_a[58]*stride_image + i_x_a[58]*n_channels); b[58]= *(image_center + i_y_b[58]*stride_image + i_x_b[58]*n_channels);
            a[59]= *(image_center + i_y_a[59]*stride_image + i_x_a[59]*n_channels); b[59]= *(image_center + i_y_b[59]*stride_image + i_x_b[59]*n_channels);
            a[60]= *(image_center + i_y_a[60]*stride_image + i_x_a[60]*n_channels); b[60]= *(image_center + i_y_b[60]*stride_image + i_x_b[60]*n_channels);
            a[61]= *(image_center + i_y_a[61]*stride_image + i_x_a[61]*n_channels); b[61]= *(image_center + i_y_b[61]*stride_image + i_x_b[61]*n_channels);
            a[62]= *(image_center + i_y_a[62]*stride_image + i_x_a[62]*n_channels); b[62]= *(image_center + i_y_b[62]*stride_image + i_x_b[62]*n_channels);
            a[63]= *(image_center + i_y_a[63]*stride_image + i_x_a[63]*n_channels); b[63]= *(image_center + i_y_b[63]*stride_image + i_x_b[63]*n_channels);
            a[64]= *(image_center + i_y_a[64]*stride_image + i_x_a[64]*n_channels); b[64]= *(image_center + i_y_b[64]*stride_image + i_x_b[64]*n_channels);
            a[65]= *(image_center + i_y_a[65]*stride_image + i_x_a[65]*n_channels); b[65]= *(image_center + i_y_b[65]*stride_image + i_x_b[65]*n_channels);
            a[66]= *(image_center + i_y_a[66]*stride_image + i_x_a[66]*n_channels); b[66]= *(image_center + i_y_b[66]*stride_image + i_x_b[66]*n_channels);
            a[67]= *(image_center + i_y_a[67]*stride_image + i_x_a[67]*n_channels); b[67]= *(image_center + i_y_b[67]*stride_image + i_x_b[67]*n_channels);
            a[68]= *(image_center + i_y_a[68]*stride_image + i_x_a[68]*n_channels); b[68]= *(image_center + i_y_b[68]*stride_image + i_x_b[68]*n_channels);
            a[69]= *(image_center + i_y_a[69]*stride_image + i_x_a[69]*n_channels); b[69]= *(image_center + i_y_b[69]*stride_image + i_x_b[69]*n_channels);
            a[70]= *(image_center + i_y_a[70]*stride_image + i_x_a[70]*n_channels); b[70]= *(image_center + i_y_b[70]*stride_image + i_x_b[70]*n_channels);
            a[71]= *(image_center + i_y_a[71]*stride_image + i_x_a[71]*n_channels); b[71]= *(image_center + i_y_b[71]*stride_image + i_x_b[71]*n_channels);
            a[72]= *(image_center + i_y_a[72]*stride_image + i_x_a[72]*n_channels); b[72]= *(image_center + i_y_b[72]*stride_image + i_x_b[72]*n_channels);
            a[73]= *(image_center + i_y_a[73]*stride_image + i_x_a[73]*n_channels); b[73]= *(image_center + i_y_b[73]*stride_image + i_x_b[73]*n_channels);
            a[74]= *(image_center + i_y_a[74]*stride_image + i_x_a[74]*n_channels); b[74]= *(image_center + i_y_b[74]*stride_image + i_x_b[74]*n_channels);
            a[75]= *(image_center + i_y_a[75]*stride_image + i_x_a[75]*n_channels); b[75]= *(image_center + i_y_b[75]*stride_image + i_x_b[75]*n_channels);
            a[76]= *(image_center + i_y_a[76]*stride_image + i_x_a[76]*n_channels); b[76]= *(image_center + i_y_b[76]*stride_image + i_x_b[76]*n_channels);
            a[77]= *(image_center + i_y_a[77]*stride_image + i_x_a[77]*n_channels); b[77]= *(image_center + i_y_b[77]*stride_image + i_x_b[77]*n_channels);
            a[78]= *(image_center + i_y_a[78]*stride_image + i_x_a[78]*n_channels); b[78]= *(image_center + i_y_b[78]*stride_image + i_x_b[78]*n_channels);
            a[79]= *(image_center + i_y_a[79]*stride_image + i_x_a[79]*n_channels); b[79]= *(image_center + i_y_b[79]*stride_image + i_x_b[79]*n_channels);
            a[80]= *(image_center + i_y_a[80]*stride_image + i_x_a[80]*n_channels); b[80]= *(image_center + i_y_b[80]*stride_image + i_x_b[80]*n_channels);
            a[81]= *(image_center + i_y_a[81]*stride_image + i_x_a[81]*n_channels); b[81]= *(image_center + i_y_b[81]*stride_image + i_x_b[81]*n_channels);
            a[82]= *(image_center + i_y_a[82]*stride_image + i_x_a[82]*n_channels); b[82]= *(image_center + i_y_b[82]*stride_image + i_x_b[82]*n_channels);
            a[83]= *(image_center + i_y_a[83]*stride_image + i_x_a[83]*n_channels); b[83]= *(image_center + i_y_b[83]*stride_image + i_x_b[83]*n_channels);
            a[84]= *(image_center + i_y_a[84]*stride_image + i_x_a[84]*n_channels); b[84]= *(image_center + i_y_b[84]*stride_image + i_x_b[84]*n_channels);
            a[85]= *(image_center + i_y_a[85]*stride_image + i_x_a[85]*n_channels); b[85]= *(image_center + i_y_b[85]*stride_image + i_x_b[85]*n_channels);
            a[86]= *(image_center + i_y_a[86]*stride_image + i_x_a[86]*n_channels); b[86]= *(image_center + i_y_b[86]*stride_image + i_x_b[86]*n_channels);
            a[87]= *(image_center + i_y_a[87]*stride_image + i_x_a[87]*n_channels); b[87]= *(image_center + i_y_b[87]*stride_image + i_x_b[87]*n_channels);
            a[88]= *(image_center + i_y_a[88]*stride_image + i_x_a[88]*n_channels); b[88]= *(image_center + i_y_b[88]*stride_image + i_x_b[88]*n_channels);
            a[89]= *(image_center + i_y_a[89]*stride_image + i_x_a[89]*n_channels); b[89]= *(image_center + i_y_b[89]*stride_image + i_x_b[89]*n_channels);
            a[90]= *(image_center + i_y_a[90]*stride_image + i_x_a[90]*n_channels); b[90]= *(image_center + i_y_b[90]*stride_image + i_x_b[90]*n_channels);
            a[91]= *(image_center + i_y_a[91]*stride_image + i_x_a[91]*n_channels); b[91]= *(image_center + i_y_b[91]*stride_image + i_x_b[91]*n_channels);
            a[92]= *(image_center + i_y_a[92]*stride_image + i_x_a[92]*n_channels); b[92]= *(image_center + i_y_b[92]*stride_image + i_x_b[92]*n_channels);
            a[93]= *(image_center + i_y_a[93]*stride_image + i_x_a[93]*n_channels); b[93]= *(image_center + i_y_b[93]*stride_image + i_x_b[93]*n_channels);
            a[94]= *(image_center + i_y_a[94]*stride_image + i_x_a[94]*n_channels); b[94]= *(image_center + i_y_b[94]*stride_image + i_x_b[94]*n_channels);
            a[95]= *(image_center + i_y_a[95]*stride_image + i_x_a[95]*n_channels); b[95]= *(image_center + i_y_b[95]*stride_image + i_x_b[95]*n_channels);
            a[96]= *(image_center + i_y_a[96]*stride_image + i_x_a[96]*n_channels); b[96]= *(image_center + i_y_b[96]*stride_image + i_x_b[96]*n_channels);
            a[97]= *(image_center + i_y_a[97]*stride_image + i_x_a[97]*n_channels); b[97]= *(image_center + i_y_b[97]*stride_image + i_x_b[97]*n_channels);
            a[98]= *(image_center + i_y_a[98]*stride_image + i_x_a[98]*n_channels); b[98]= *(image_center + i_y_b[98]*stride_image + i_x_b[98]*n_channels);
            a[99]= *(image_center + i_y_a[99]*stride_image + i_x_a[99]*n_channels); b[99]= *(image_center + i_y_b[99]*stride_image + i_x_b[99]*n_channels);
            a[100]= *(image_center + i_y_a[100]*stride_image + i_x_a[100]*n_channels); b[100]= *(image_center + i_y_b[100]*stride_image + i_x_b[100]*n_channels);
            a[101]= *(image_center + i_y_a[101]*stride_image + i_x_a[101]*n_channels); b[101]= *(image_center + i_y_b[101]*stride_image + i_x_b[101]*n_channels);
            a[102]= *(image_center + i_y_a[102]*stride_image + i_x_a[102]*n_channels); b[102]= *(image_center + i_y_b[102]*stride_image + i_x_b[102]*n_channels);
            a[103]= *(image_center + i_y_a[103]*stride_image + i_x_a[103]*n_channels); b[103]= *(image_center + i_y_b[103]*stride_image + i_x_b[103]*n_channels);
            a[104]= *(image_center + i_y_a[104]*stride_image + i_x_a[104]*n_channels); b[104]= *(image_center + i_y_b[104]*stride_image + i_x_b[104]*n_channels);
            a[105]= *(image_center + i_y_a[105]*stride_image + i_x_a[105]*n_channels); b[105]= *(image_center + i_y_b[105]*stride_image + i_x_b[105]*n_channels);
            a[106]= *(image_center + i_y_a[106]*stride_image + i_x_a[106]*n_channels); b[106]= *(image_center + i_y_b[106]*stride_image + i_x_b[106]*n_channels);
            a[107]= *(image_center + i_y_a[107]*stride_image + i_x_a[107]*n_channels); b[107]= *(image_center + i_y_b[107]*stride_image + i_x_b[107]*n_channels);
            a[108]= *(image_center + i_y_a[108]*stride_image + i_x_a[108]*n_channels); b[108]= *(image_center + i_y_b[108]*stride_image + i_x_b[108]*n_channels);
            a[109]= *(image_center + i_y_a[109]*stride_image + i_x_a[109]*n_channels); b[109]= *(image_center + i_y_b[109]*stride_image + i_x_b[109]*n_channels);
            a[110]= *(image_center + i_y_a[110]*stride_image + i_x_a[110]*n_channels); b[110]= *(image_center + i_y_b[110]*stride_image + i_x_b[110]*n_channels);
            a[111]= *(image_center + i_y_a[111]*stride_image + i_x_a[111]*n_channels); b[111]= *(image_center + i_y_b[111]*stride_image + i_x_b[111]*n_channels);
            a[112]= *(image_center + i_y_a[112]*stride_image + i_x_a[112]*n_channels); b[112]= *(image_center + i_y_b[112]*stride_image + i_x_b[112]*n_channels);
            a[113]= *(image_center + i_y_a[113]*stride_image + i_x_a[113]*n_channels); b[113]= *(image_center + i_y_b[113]*stride_image + i_x_b[113]*n_channels);
            a[114]= *(image_center + i_y_a[114]*stride_image + i_x_a[114]*n_channels); b[114]= *(image_center + i_y_b[114]*stride_image + i_x_b[114]*n_channels);
            a[115]= *(image_center + i_y_a[115]*stride_image + i_x_a[115]*n_channels); b[115]= *(image_center + i_y_b[115]*stride_image + i_x_b[115]*n_channels);
            a[116]= *(image_center + i_y_a[116]*stride_image + i_x_a[116]*n_channels); b[116]= *(image_center + i_y_b[116]*stride_image + i_x_b[116]*n_channels);
            a[117]= *(image_center + i_y_a[117]*stride_image + i_x_a[117]*n_channels); b[117]= *(image_center + i_y_b[117]*stride_image + i_x_b[117]*n_channels);
            a[118]= *(image_center + i_y_a[118]*stride_image + i_x_a[118]*n_channels); b[118]= *(image_center + i_y_b[118]*stride_image + i_x_b[118]*n_channels);
            a[119]= *(image_center + i_y_a[119]*stride_image + i_x_a[119]*n_channels); b[119]= *(image_center + i_y_b[119]*stride_image + i_x_b[119]*n_channels);
            a[120]= *(image_center + i_y_a[120]*stride_image + i_x_a[120]*n_channels); b[120]= *(image_center + i_y_b[120]*stride_image + i_x_b[120]*n_channels);
            a[121]= *(image_center + i_y_a[121]*stride_image + i_x_a[121]*n_channels); b[121]= *(image_center + i_y_b[121]*stride_image + i_x_b[121]*n_channels);
            a[122]= *(image_center + i_y_a[122]*stride_image + i_x_a[122]*n_channels); b[122]= *(image_center + i_y_b[122]*stride_image + i_x_b[122]*n_channels);
            a[123]= *(image_center + i_y_a[123]*stride_image + i_x_a[123]*n_channels); b[123]= *(image_center + i_y_b[123]*stride_image + i_x_b[123]*n_channels);
            a[124]= *(image_center + i_y_a[124]*stride_image + i_x_a[124]*n_channels); b[124]= *(image_center + i_y_b[124]*stride_image + i_x_b[124]*n_channels);
            a[125]= *(image_center + i_y_a[125]*stride_image + i_x_a[125]*n_channels); b[125]= *(image_center + i_y_b[125]*stride_image + i_x_b[125]*n_channels);
            a[126]= *(image_center + i_y_a[126]*stride_image + i_x_a[126]*n_channels); b[126]= *(image_center + i_y_b[126]*stride_image + i_x_b[126]*n_channels);
            a[127]= *(image_center + i_y_a[127]*stride_image + i_x_a[127]*n_channels); b[127]= *(image_center + i_y_b[127]*stride_image + i_x_b[127]*n_channels);
            a[128]= *(image_center + i_y_a[128]*stride_image + i_x_a[128]*n_channels); b[128]= *(image_center + i_y_b[128]*stride_image + i_x_b[128]*n_channels);
            a[129]= *(image_center + i_y_a[129]*stride_image + i_x_a[129]*n_channels); b[129]= *(image_center + i_y_b[129]*stride_image + i_x_b[129]*n_channels);
            a[130]= *(image_center + i_y_a[130]*stride_image + i_x_a[130]*n_channels); b[130]= *(image_center + i_y_b[130]*stride_image + i_x_b[130]*n_channels);
            a[131]= *(image_center + i_y_a[131]*stride_image + i_x_a[131]*n_channels); b[131]= *(image_center + i_y_b[131]*stride_image + i_x_b[131]*n_channels);
            a[132]= *(image_center + i_y_a[132]*stride_image + i_x_a[132]*n_channels); b[132]= *(image_center + i_y_b[132]*stride_image + i_x_b[132]*n_channels);
            a[133]= *(image_center + i_y_a[133]*stride_image + i_x_a[133]*n_channels); b[133]= *(image_center + i_y_b[133]*stride_image + i_x_b[133]*n_channels);
            a[134]= *(image_center + i_y_a[134]*stride_image + i_x_a[134]*n_channels); b[134]= *(image_center + i_y_b[134]*stride_image + i_x_b[134]*n_channels);
            a[135]= *(image_center + i_y_a[135]*stride_image + i_x_a[135]*n_channels); b[135]= *(image_center + i_y_b[135]*stride_image + i_x_b[135]*n_channels);
            a[136]= *(image_center + i_y_a[136]*stride_image + i_x_a[136]*n_channels); b[136]= *(image_center + i_y_b[136]*stride_image + i_x_b[136]*n_channels);
            a[137]= *(image_center + i_y_a[137]*stride_image + i_x_a[137]*n_channels); b[137]= *(image_center + i_y_b[137]*stride_image + i_x_b[137]*n_channels);
            a[138]= *(image_center + i_y_a[138]*stride_image + i_x_a[138]*n_channels); b[138]= *(image_center + i_y_b[138]*stride_image + i_x_b[138]*n_channels);
            a[139]= *(image_center + i_y_a[139]*stride_image + i_x_a[139]*n_channels); b[139]= *(image_center + i_y_b[139]*stride_image + i_x_b[139]*n_channels);
            a[140]= *(image_center + i_y_a[140]*stride_image + i_x_a[140]*n_channels); b[140]= *(image_center + i_y_b[140]*stride_image + i_x_b[140]*n_channels);
            a[141]= *(image_center + i_y_a[141]*stride_image + i_x_a[141]*n_channels); b[141]= *(image_center + i_y_b[141]*stride_image + i_x_b[141]*n_channels);
            a[142]= *(image_center + i_y_a[142]*stride_image + i_x_a[142]*n_channels); b[142]= *(image_center + i_y_b[142]*stride_image + i_x_b[142]*n_channels);
            a[143]= *(image_center + i_y_a[143]*stride_image + i_x_a[143]*n_channels); b[143]= *(image_center + i_y_b[143]*stride_image + i_x_b[143]*n_channels);
            a[144]= *(image_center + i_y_a[144]*stride_image + i_x_a[144]*n_channels); b[144]= *(image_center + i_y_b[144]*stride_image + i_x_b[144]*n_channels);
            a[145]= *(image_center + i_y_a[145]*stride_image + i_x_a[145]*n_channels); b[145]= *(image_center + i_y_b[145]*stride_image + i_x_b[145]*n_channels);
            a[146]= *(image_center + i_y_a[146]*stride_image + i_x_a[146]*n_channels); b[146]= *(image_center + i_y_b[146]*stride_image + i_x_b[146]*n_channels);
            a[147]= *(image_center + i_y_a[147]*stride_image + i_x_a[147]*n_channels); b[147]= *(image_center + i_y_b[147]*stride_image + i_x_b[147]*n_channels);
            a[148]= *(image_center + i_y_a[148]*stride_image + i_x_a[148]*n_channels); b[148]= *(image_center + i_y_b[148]*stride_image + i_x_b[148]*n_channels);
            a[149]= *(image_center + i_y_a[149]*stride_image + i_x_a[149]*n_channels); b[149]= *(image_center + i_y_b[149]*stride_image + i_x_b[149]*n_channels);
            a[150]= *(image_center + i_y_a[150]*stride_image + i_x_a[150]*n_channels); b[150]= *(image_center + i_y_b[150]*stride_image + i_x_b[150]*n_channels);
            a[151]= *(image_center + i_y_a[151]*stride_image + i_x_a[151]*n_channels); b[151]= *(image_center + i_y_b[151]*stride_image + i_x_b[151]*n_channels);
            a[152]= *(image_center + i_y_a[152]*stride_image + i_x_a[152]*n_channels); b[152]= *(image_center + i_y_b[152]*stride_image + i_x_b[152]*n_channels);
            a[153]= *(image_center + i_y_a[153]*stride_image + i_x_a[153]*n_channels); b[153]= *(image_center + i_y_b[153]*stride_image + i_x_b[153]*n_channels);
            a[154]= *(image_center + i_y_a[154]*stride_image + i_x_a[154]*n_channels); b[154]= *(image_center + i_y_b[154]*stride_image + i_x_b[154]*n_channels);
            a[155]= *(image_center + i_y_a[155]*stride_image + i_x_a[155]*n_channels); b[155]= *(image_center + i_y_b[155]*stride_image + i_x_b[155]*n_channels);
            a[156]= *(image_center + i_y_a[156]*stride_image + i_x_a[156]*n_channels); b[156]= *(image_center + i_y_b[156]*stride_image + i_x_b[156]*n_channels);
            a[157]= *(image_center + i_y_a[157]*stride_image + i_x_a[157]*n_channels); b[157]= *(image_center + i_y_b[157]*stride_image + i_x_b[157]*n_channels);
            a[158]= *(image_center + i_y_a[158]*stride_image + i_x_a[158]*n_channels); b[158]= *(image_center + i_y_b[158]*stride_image + i_x_b[158]*n_channels);
            a[159]= *(image_center + i_y_a[159]*stride_image + i_x_a[159]*n_channels); b[159]= *(image_center + i_y_b[159]*stride_image + i_x_b[159]*n_channels);
            a[160]= *(image_center + i_y_a[160]*stride_image + i_x_a[160]*n_channels); b[160]= *(image_center + i_y_b[160]*stride_image + i_x_b[160]*n_channels);
            a[161]= *(image_center + i_y_a[161]*stride_image + i_x_a[161]*n_channels); b[161]= *(image_center + i_y_b[161]*stride_image + i_x_b[161]*n_channels);
            a[162]= *(image_center + i_y_a[162]*stride_image + i_x_a[162]*n_channels); b[162]= *(image_center + i_y_b[162]*stride_image + i_x_b[162]*n_channels);
            a[163]= *(image_center + i_y_a[163]*stride_image + i_x_a[163]*n_channels); b[163]= *(image_center + i_y_b[163]*stride_image + i_x_b[163]*n_channels);
            a[164]= *(image_center + i_y_a[164]*stride_image + i_x_a[164]*n_channels); b[164]= *(image_center + i_y_b[164]*stride_image + i_x_b[164]*n_channels);
            a[165]= *(image_center + i_y_a[165]*stride_image + i_x_a[165]*n_channels); b[165]= *(image_center + i_y_b[165]*stride_image + i_x_b[165]*n_channels);
            a[166]= *(image_center + i_y_a[166]*stride_image + i_x_a[166]*n_channels); b[166]= *(image_center + i_y_b[166]*stride_image + i_x_b[166]*n_channels);
            a[167]= *(image_center + i_y_a[167]*stride_image + i_x_a[167]*n_channels); b[167]= *(image_center + i_y_b[167]*stride_image + i_x_b[167]*n_channels);
            a[168]= *(image_center + i_y_a[168]*stride_image + i_x_a[168]*n_channels); b[168]= *(image_center + i_y_b[168]*stride_image + i_x_b[168]*n_channels);
            a[169]= *(image_center + i_y_a[169]*stride_image + i_x_a[169]*n_channels); b[169]= *(image_center + i_y_b[169]*stride_image + i_x_b[169]*n_channels);
            a[170]= *(image_center + i_y_a[170]*stride_image + i_x_a[170]*n_channels); b[170]= *(image_center + i_y_b[170]*stride_image + i_x_b[170]*n_channels);
            a[171]= *(image_center + i_y_a[171]*stride_image + i_x_a[171]*n_channels); b[171]= *(image_center + i_y_b[171]*stride_image + i_x_b[171]*n_channels);
            a[172]= *(image_center + i_y_a[172]*stride_image + i_x_a[172]*n_channels); b[172]= *(image_center + i_y_b[172]*stride_image + i_x_b[172]*n_channels);
            a[173]= *(image_center + i_y_a[173]*stride_image + i_x_a[173]*n_channels); b[173]= *(image_center + i_y_b[173]*stride_image + i_x_b[173]*n_channels);
            a[174]= *(image_center + i_y_a[174]*stride_image + i_x_a[174]*n_channels); b[174]= *(image_center + i_y_b[174]*stride_image + i_x_b[174]*n_channels);
            a[175]= *(image_center + i_y_a[175]*stride_image + i_x_a[175]*n_channels); b[175]= *(image_center + i_y_b[175]*stride_image + i_x_b[175]*n_channels);
            a[176]= *(image_center + i_y_a[176]*stride_image + i_x_a[176]*n_channels); b[176]= *(image_center + i_y_b[176]*stride_image + i_x_b[176]*n_channels);
            a[177]= *(image_center + i_y_a[177]*stride_image + i_x_a[177]*n_channels); b[177]= *(image_center + i_y_b[177]*stride_image + i_x_b[177]*n_channels);
            a[178]= *(image_center + i_y_a[178]*stride_image + i_x_a[178]*n_channels); b[178]= *(image_center + i_y_b[178]*stride_image + i_x_b[178]*n_channels);
            a[179]= *(image_center + i_y_a[179]*stride_image + i_x_a[179]*n_channels); b[179]= *(image_center + i_y_b[179]*stride_image + i_x_b[179]*n_channels);
            a[180]= *(image_center + i_y_a[180]*stride_image + i_x_a[180]*n_channels); b[180]= *(image_center + i_y_b[180]*stride_image + i_x_b[180]*n_channels);
            a[181]= *(image_center + i_y_a[181]*stride_image + i_x_a[181]*n_channels); b[181]= *(image_center + i_y_b[181]*stride_image + i_x_b[181]*n_channels);
            a[182]= *(image_center + i_y_a[182]*stride_image + i_x_a[182]*n_channels); b[182]= *(image_center + i_y_b[182]*stride_image + i_x_b[182]*n_channels);
            a[183]= *(image_center + i_y_a[183]*stride_image + i_x_a[183]*n_channels); b[183]= *(image_center + i_y_b[183]*stride_image + i_x_b[183]*n_channels);
            a[184]= *(image_center + i_y_a[184]*stride_image + i_x_a[184]*n_channels); b[184]= *(image_center + i_y_b[184]*stride_image + i_x_b[184]*n_channels);
            a[185]= *(image_center + i_y_a[185]*stride_image + i_x_a[185]*n_channels); b[185]= *(image_center + i_y_b[185]*stride_image + i_x_b[185]*n_channels);
            a[186]= *(image_center + i_y_a[186]*stride_image + i_x_a[186]*n_channels); b[186]= *(image_center + i_y_b[186]*stride_image + i_x_b[186]*n_channels);
            a[187]= *(image_center + i_y_a[187]*stride_image + i_x_a[187]*n_channels); b[187]= *(image_center + i_y_b[187]*stride_image + i_x_b[187]*n_channels);
            a[188]= *(image_center + i_y_a[188]*stride_image + i_x_a[188]*n_channels); b[188]= *(image_center + i_y_b[188]*stride_image + i_x_b[188]*n_channels);
            a[189]= *(image_center + i_y_a[189]*stride_image + i_x_a[189]*n_channels); b[189]= *(image_center + i_y_b[189]*stride_image + i_x_b[189]*n_channels);
            a[190]= *(image_center + i_y_a[190]*stride_image + i_x_a[190]*n_channels); b[190]= *(image_center + i_y_b[190]*stride_image + i_x_b[190]*n_channels);
            a[191]= *(image_center + i_y_a[191]*stride_image + i_x_a[191]*n_channels); b[191]= *(image_center + i_y_b[191]*stride_image + i_x_b[191]*n_channels);
            a[192]= *(image_center + i_y_a[192]*stride_image + i_x_a[192]*n_channels); b[192]= *(image_center + i_y_b[192]*stride_image + i_x_b[192]*n_channels);
            a[193]= *(image_center + i_y_a[193]*stride_image + i_x_a[193]*n_channels); b[193]= *(image_center + i_y_b[193]*stride_image + i_x_b[193]*n_channels);
            a[194]= *(image_center + i_y_a[194]*stride_image + i_x_a[194]*n_channels); b[194]= *(image_center + i_y_b[194]*stride_image + i_x_b[194]*n_channels);
            a[195]= *(image_center + i_y_a[195]*stride_image + i_x_a[195]*n_channels); b[195]= *(image_center + i_y_b[195]*stride_image + i_x_b[195]*n_channels);
            a[196]= *(image_center + i_y_a[196]*stride_image + i_x_a[196]*n_channels); b[196]= *(image_center + i_y_b[196]*stride_image + i_x_b[196]*n_channels);
            a[197]= *(image_center + i_y_a[197]*stride_image + i_x_a[197]*n_channels); b[197]= *(image_center + i_y_b[197]*stride_image + i_x_b[197]*n_channels);
            a[198]= *(image_center + i_y_a[198]*stride_image + i_x_a[198]*n_channels); b[198]= *(image_center + i_y_b[198]*stride_image + i_x_b[198]*n_channels);
            a[199]= *(image_center + i_y_a[199]*stride_image + i_x_a[199]*n_channels); b[199]= *(image_center + i_y_b[199]*stride_image + i_x_b[199]*n_channels);
            a[200]= *(image_center + i_y_a[200]*stride_image + i_x_a[200]*n_channels); b[200]= *(image_center + i_y_b[200]*stride_image + i_x_b[200]*n_channels);
            a[201]= *(image_center + i_y_a[201]*stride_image + i_x_a[201]*n_channels); b[201]= *(image_center + i_y_b[201]*stride_image + i_x_b[201]*n_channels);
            a[202]= *(image_center + i_y_a[202]*stride_image + i_x_a[202]*n_channels); b[202]= *(image_center + i_y_b[202]*stride_image + i_x_b[202]*n_channels);
            a[203]= *(image_center + i_y_a[203]*stride_image + i_x_a[203]*n_channels); b[203]= *(image_center + i_y_b[203]*stride_image + i_x_b[203]*n_channels);
            a[204]= *(image_center + i_y_a[204]*stride_image + i_x_a[204]*n_channels); b[204]= *(image_center + i_y_b[204]*stride_image + i_x_b[204]*n_channels);
            a[205]= *(image_center + i_y_a[205]*stride_image + i_x_a[205]*n_channels); b[205]= *(image_center + i_y_b[205]*stride_image + i_x_b[205]*n_channels);
            a[206]= *(image_center + i_y_a[206]*stride_image + i_x_a[206]*n_channels); b[206]= *(image_center + i_y_b[206]*stride_image + i_x_b[206]*n_channels);
            a[207]= *(image_center + i_y_a[207]*stride_image + i_x_a[207]*n_channels); b[207]= *(image_center + i_y_b[207]*stride_image + i_x_b[207]*n_channels);
            a[208]= *(image_center + i_y_a[208]*stride_image + i_x_a[208]*n_channels); b[208]= *(image_center + i_y_b[208]*stride_image + i_x_b[208]*n_channels);
            a[209]= *(image_center + i_y_a[209]*stride_image + i_x_a[209]*n_channels); b[209]= *(image_center + i_y_b[209]*stride_image + i_x_b[209]*n_channels);
            a[210]= *(image_center + i_y_a[210]*stride_image + i_x_a[210]*n_channels); b[210]= *(image_center + i_y_b[210]*stride_image + i_x_b[210]*n_channels);
            a[211]= *(image_center + i_y_a[211]*stride_image + i_x_a[211]*n_channels); b[211]= *(image_center + i_y_b[211]*stride_image + i_x_b[211]*n_channels);
            a[212]= *(image_center + i_y_a[212]*stride_image + i_x_a[212]*n_channels); b[212]= *(image_center + i_y_b[212]*stride_image + i_x_b[212]*n_channels);
            a[213]= *(image_center + i_y_a[213]*stride_image + i_x_a[213]*n_channels); b[213]= *(image_center + i_y_b[213]*stride_image + i_x_b[213]*n_channels);
            a[214]= *(image_center + i_y_a[214]*stride_image + i_x_a[214]*n_channels); b[214]= *(image_center + i_y_b[214]*stride_image + i_x_b[214]*n_channels);
            a[215]= *(image_center + i_y_a[215]*stride_image + i_x_a[215]*n_channels); b[215]= *(image_center + i_y_b[215]*stride_image + i_x_b[215]*n_channels);
            a[216]= *(image_center + i_y_a[216]*stride_image + i_x_a[216]*n_channels); b[216]= *(image_center + i_y_b[216]*stride_image + i_x_b[216]*n_channels);
            a[217]= *(image_center + i_y_a[217]*stride_image + i_x_a[217]*n_channels); b[217]= *(image_center + i_y_b[217]*stride_image + i_x_b[217]*n_channels);
            a[218]= *(image_center + i_y_a[218]*stride_image + i_x_a[218]*n_channels); b[218]= *(image_center + i_y_b[218]*stride_image + i_x_b[218]*n_channels);
            a[219]= *(image_center + i_y_a[219]*stride_image + i_x_a[219]*n_channels); b[219]= *(image_center + i_y_b[219]*stride_image + i_x_b[219]*n_channels);
            a[220]= *(image_center + i_y_a[220]*stride_image + i_x_a[220]*n_channels); b[220]= *(image_center + i_y_b[220]*stride_image + i_x_b[220]*n_channels);
            a[221]= *(image_center + i_y_a[221]*stride_image + i_x_a[221]*n_channels); b[221]= *(image_center + i_y_b[221]*stride_image + i_x_b[221]*n_channels);
            a[222]= *(image_center + i_y_a[222]*stride_image + i_x_a[222]*n_channels); b[222]= *(image_center + i_y_b[222]*stride_image + i_x_b[222]*n_channels);
            a[223]= *(image_center + i_y_a[223]*stride_image + i_x_a[223]*n_channels); b[223]= *(image_center + i_y_b[223]*stride_image + i_x_b[223]*n_channels);
            a[224]= *(image_center + i_y_a[224]*stride_image + i_x_a[224]*n_channels); b[224]= *(image_center + i_y_b[224]*stride_image + i_x_b[224]*n_channels);
            a[225]= *(image_center + i_y_a[225]*stride_image + i_x_a[225]*n_channels); b[225]= *(image_center + i_y_b[225]*stride_image + i_x_b[225]*n_channels);
            a[226]= *(image_center + i_y_a[226]*stride_image + i_x_a[226]*n_channels); b[226]= *(image_center + i_y_b[226]*stride_image + i_x_b[226]*n_channels);
            a[227]= *(image_center + i_y_a[227]*stride_image + i_x_a[227]*n_channels); b[227]= *(image_center + i_y_b[227]*stride_image + i_x_b[227]*n_channels);
            a[228]= *(image_center + i_y_a[228]*stride_image + i_x_a[228]*n_channels); b[228]= *(image_center + i_y_b[228]*stride_image + i_x_b[228]*n_channels);
            a[229]= *(image_center + i_y_a[229]*stride_image + i_x_a[229]*n_channels); b[229]= *(image_center + i_y_b[229]*stride_image + i_x_b[229]*n_channels);
            a[230]= *(image_center + i_y_a[230]*stride_image + i_x_a[230]*n_channels); b[230]= *(image_center + i_y_b[230]*stride_image + i_x_b[230]*n_channels);
            a[231]= *(image_center + i_y_a[231]*stride_image + i_x_a[231]*n_channels); b[231]= *(image_center + i_y_b[231]*stride_image + i_x_b[231]*n_channels);
            a[232]= *(image_center + i_y_a[232]*stride_image + i_x_a[232]*n_channels); b[232]= *(image_center + i_y_b[232]*stride_image + i_x_b[232]*n_channels);
            a[233]= *(image_center + i_y_a[233]*stride_image + i_x_a[233]*n_channels); b[233]= *(image_center + i_y_b[233]*stride_image + i_x_b[233]*n_channels);
            a[234]= *(image_center + i_y_a[234]*stride_image + i_x_a[234]*n_channels); b[234]= *(image_center + i_y_b[234]*stride_image + i_x_b[234]*n_channels);
            a[235]= *(image_center + i_y_a[235]*stride_image + i_x_a[235]*n_channels); b[235]= *(image_center + i_y_b[235]*stride_image + i_x_b[235]*n_channels);
            a[236]= *(image_center + i_y_a[236]*stride_image + i_x_a[236]*n_channels); b[236]= *(image_center + i_y_b[236]*stride_image + i_x_b[236]*n_channels);
            a[237]= *(image_center + i_y_a[237]*stride_image + i_x_a[237]*n_channels); b[237]= *(image_center + i_y_b[237]*stride_image + i_x_b[237]*n_channels);
            a[238]= *(image_center + i_y_a[238]*stride_image + i_x_a[238]*n_channels); b[238]= *(image_center + i_y_b[238]*stride_image + i_x_b[238]*n_channels);
            a[239]= *(image_center + i_y_a[239]*stride_image + i_x_a[239]*n_channels); b[239]= *(image_center + i_y_b[239]*stride_image + i_x_b[239]*n_channels);
            a[240]= *(image_center + i_y_a[240]*stride_image + i_x_a[240]*n_channels); b[240]= *(image_center + i_y_b[240]*stride_image + i_x_b[240]*n_channels);
            a[241]= *(image_center + i_y_a[241]*stride_image + i_x_a[241]*n_channels); b[241]= *(image_center + i_y_b[241]*stride_image + i_x_b[241]*n_channels);
            a[242]= *(image_center + i_y_a[242]*stride_image + i_x_a[242]*n_channels); b[242]= *(image_center + i_y_b[242]*stride_image + i_x_b[242]*n_channels);
            a[243]= *(image_center + i_y_a[243]*stride_image + i_x_a[243]*n_channels); b[243]= *(image_center + i_y_b[243]*stride_image + i_x_b[243]*n_channels);
            a[244]= *(image_center + i_y_a[244]*stride_image + i_x_a[244]*n_channels); b[244]= *(image_center + i_y_b[244]*stride_image + i_x_b[244]*n_channels);
            a[245]= *(image_center + i_y_a[245]*stride_image + i_x_a[245]*n_channels); b[245]= *(image_center + i_y_b[245]*stride_image + i_x_b[245]*n_channels);
            a[246]= *(image_center + i_y_a[246]*stride_image + i_x_a[246]*n_channels); b[246]= *(image_center + i_y_b[246]*stride_image + i_x_b[246]*n_channels);
            a[247]= *(image_center + i_y_a[247]*stride_image + i_x_a[247]*n_channels); b[247]= *(image_center + i_y_b[247]*stride_image + i_x_b[247]*n_channels);
            a[248]= *(image_center + i_y_a[248]*stride_image + i_x_a[248]*n_channels); b[248]= *(image_center + i_y_b[248]*stride_image + i_x_b[248]*n_channels);
            a[249]= *(image_center + i_y_a[249]*stride_image + i_x_a[249]*n_channels); b[249]= *(image_center + i_y_b[249]*stride_image + i_x_b[249]*n_channels);
            a[250]= *(image_center + i_y_a[250]*stride_image + i_x_a[250]*n_channels); b[250]= *(image_center + i_y_b[250]*stride_image + i_x_b[250]*n_channels);
            a[251]= *(image_center + i_y_a[251]*stride_image + i_x_a[251]*n_channels); b[251]= *(image_center + i_y_b[251]*stride_image + i_x_b[251]*n_channels);
            a[252]= *(image_center + i_y_a[252]*stride_image + i_x_a[252]*n_channels); b[252]= *(image_center + i_y_b[252]*stride_image + i_x_b[252]*n_channels);
            a[253]= *(image_center + i_y_a[253]*stride_image + i_x_a[253]*n_channels); b[253]= *(image_center + i_y_b[253]*stride_image + i_x_b[253]*n_channels);
            a[254]= *(image_center + i_y_a[254]*stride_image + i_x_a[254]*n_channels); b[254]= *(image_center + i_y_b[254]*stride_image + i_x_b[254]*n_channels);
            a[255]= *(image_center + i_y_a[255]*stride_image + i_x_a[255]*n_channels); b[255]= *(image_center + i_y_b[255]*stride_image + i_x_b[255]*n_channels);
            
            
            unsigned char f[32] __attribute__((aligned(32)));
            f[0] |= (unsigned char) (a[0*8 + 0] > b[0*8 + 0]) << 0;
            f[0] |= (unsigned char) (a[0*8 + 1] > b[0*8 + 1]) << 1;
            f[0] |= (unsigned char) (a[0*8 + 2] > b[0*8 + 2]) << 2;
            f[0] |= (unsigned char) (a[0*8 + 3] > b[0*8 + 3]) << 3;
            f[0] |= (unsigned char) (a[0*8 + 4] > b[0*8 + 4]) << 4;
            f[0] |= (unsigned char) (a[0*8 + 5] > b[0*8 + 5]) << 5;
            f[0] |= (unsigned char) (a[0*8 + 6] > b[0*8 + 6]) << 6;
            f[0] |= (unsigned char) (a[0*8 + 7] > b[0*8 + 7]) << 7;
            
            f[1] |= (unsigned char) (a[1*8 + 0] > b[1*8 + 0]) << 0;
            f[1] |= (unsigned char) (a[1*8 + 1] > b[1*8 + 1]) << 1;
            f[1] |= (unsigned char) (a[1*8 + 2] > b[1*8 + 2]) << 2;
            f[1] |= (unsigned char) (a[1*8 + 3] > b[1*8 + 3]) << 3;
            f[1] |= (unsigned char) (a[1*8 + 4] > b[1*8 + 4]) << 4;
            f[1] |= (unsigned char) (a[1*8 + 5] > b[1*8 + 5]) << 5;
            f[1] |= (unsigned char) (a[1*8 + 6] > b[1*8 + 6]) << 6;
            f[1] |= (unsigned char) (a[1*8 + 7] > b[1*8 + 7]) << 7;
            
            f[2] |= (unsigned char) (a[2*8 + 0] > b[2*8 + 0]) << 0;
            f[2] |= (unsigned char) (a[2*8 + 1] > b[2*8 + 1]) << 1;
            f[2] |= (unsigned char) (a[2*8 + 2] > b[2*8 + 2]) << 2;
            f[2] |= (unsigned char) (a[2*8 + 3] > b[2*8 + 3]) << 3;
            f[2] |= (unsigned char) (a[2*8 + 4] > b[2*8 + 4]) << 4;
            f[2] |= (unsigned char) (a[2*8 + 5] > b[2*8 + 5]) << 5;
            f[2] |= (unsigned char) (a[2*8 + 6] > b[2*8 + 6]) << 6;
            f[2] |= (unsigned char) (a[2*8 + 7] > b[2*8 + 7]) << 7;
            
            f[3] |= (unsigned char) (a[3*8 + 0] > b[3*8 + 0]) << 0;
            f[3] |= (unsigned char) (a[3*8 + 1] > b[3*8 + 1]) << 1;
            f[3] |= (unsigned char) (a[3*8 + 2] > b[3*8 + 2]) << 2;
            f[3] |= (unsigned char) (a[3*8 + 3] > b[3*8 + 3]) << 3;
            f[3] |= (unsigned char) (a[3*8 + 4] > b[3*8 + 4]) << 4;
            f[3] |= (unsigned char) (a[3*8 + 5] > b[3*8 + 5]) << 5;
            f[3] |= (unsigned char) (a[3*8 + 6] > b[3*8 + 6]) << 6;
            f[3] |= (unsigned char) (a[3*8 + 7] > b[3*8 + 7]) << 7;
            
            f[4] |= (unsigned char) (a[4*8 + 0] > b[4*8 + 0]) << 0;
            f[4] |= (unsigned char) (a[4*8 + 1] > b[4*8 + 1]) << 1;
            f[4] |= (unsigned char) (a[4*8 + 2] > b[4*8 + 2]) << 2;
            f[4] |= (unsigned char) (a[4*8 + 3] > b[4*8 + 3]) << 3;
            f[4] |= (unsigned char) (a[4*8 + 4] > b[4*8 + 4]) << 4;
            f[4] |= (unsigned char) (a[4*8 + 5] > b[4*8 + 5]) << 5;
            f[4] |= (unsigned char) (a[4*8 + 6] > b[4*8 + 6]) << 6;
            f[4] |= (unsigned char) (a[4*8 + 7] > b[4*8 + 7]) << 7;
            
            f[5] |= (unsigned char) (a[5*8 + 0] > b[5*8 + 0]) << 0;
            f[5] |= (unsigned char) (a[5*8 + 1] > b[5*8 + 1]) << 1;
            f[5] |= (unsigned char) (a[5*8 + 2] > b[5*8 + 2]) << 2;
            f[5] |= (unsigned char) (a[5*8 + 3] > b[5*8 + 3]) << 3;
            f[5] |= (unsigned char) (a[5*8 + 4] > b[5*8 + 4]) << 4;
            f[5] |= (unsigned char) (a[5*8 + 5] > b[5*8 + 5]) << 5;
            f[5] |= (unsigned char) (a[5*8 + 6] > b[5*8 + 6]) << 6;
            f[5] |= (unsigned char) (a[5*8 + 7] > b[5*8 + 7]) << 7;
            
            f[6] |= (unsigned char) (a[6*8 + 0] > b[6*8 + 0]) << 0;
            f[6] |= (unsigned char) (a[6*8 + 1] > b[6*8 + 1]) << 1;
            f[6] |= (unsigned char) (a[6*8 + 2] > b[6*8 + 2]) << 2;
            f[6] |= (unsigned char) (a[6*8 + 3] > b[6*8 + 3]) << 3;
            f[6] |= (unsigned char) (a[6*8 + 4] > b[6*8 + 4]) << 4;
            f[6] |= (unsigned char) (a[6*8 + 5] > b[6*8 + 5]) << 5;
            f[6] |= (unsigned char) (a[6*8 + 6] > b[6*8 + 6]) << 6;
            f[6] |= (unsigned char) (a[6*8 + 7] > b[6*8 + 7]) << 7;
            
            f[7] |= (unsigned char) (a[7*8 + 0] > b[7*8 + 0]) << 0;
            f[7] |= (unsigned char) (a[7*8 + 1] > b[7*8 + 1]) << 1;
            f[7] |= (unsigned char) (a[7*8 + 2] > b[7*8 + 2]) << 2;
            f[7] |= (unsigned char) (a[7*8 + 3] > b[7*8 + 3]) << 3;
            f[7] |= (unsigned char) (a[7*8 + 4] > b[7*8 + 4]) << 4;
            f[7] |= (unsigned char) (a[7*8 + 5] > b[7*8 + 5]) << 5;
            f[7] |= (unsigned char) (a[7*8 + 6] > b[7*8 + 6]) << 6;
            f[7] |= (unsigned char) (a[7*8 + 7] > b[7*8 + 7]) << 7;
            
            f[8] |= (unsigned char) (a[8*8 + 0] > b[8*8 + 0]) << 0;
            f[8] |= (unsigned char) (a[8*8 + 1] > b[8*8 + 1]) << 1;
            f[8] |= (unsigned char) (a[8*8 + 2] > b[8*8 + 2]) << 2;
            f[8] |= (unsigned char) (a[8*8 + 3] > b[8*8 + 3]) << 3;
            f[8] |= (unsigned char) (a[8*8 + 4] > b[8*8 + 4]) << 4;
            f[8] |= (unsigned char) (a[8*8 + 5] > b[8*8 + 5]) << 5;
            f[8] |= (unsigned char) (a[8*8 + 6] > b[8*8 + 6]) << 6;
            f[8] |= (unsigned char) (a[8*8 + 7] > b[8*8 + 7]) << 7;
            
            f[9] |= (unsigned char) (a[9*8 + 0] > b[9*8 + 0]) << 0;
            f[9] |= (unsigned char) (a[9*8 + 1] > b[9*8 + 1]) << 1;
            f[9] |= (unsigned char) (a[9*8 + 2] > b[9*8 + 2]) << 2;
            f[9] |= (unsigned char) (a[9*8 + 3] > b[9*8 + 3]) << 3;
            f[9] |= (unsigned char) (a[9*8 + 4] > b[9*8 + 4]) << 4;
            f[9] |= (unsigned char) (a[9*8 + 5] > b[9*8 + 5]) << 5;
            f[9] |= (unsigned char) (a[9*8 + 6] > b[9*8 + 6]) << 6;
            f[9] |= (unsigned char) (a[9*8 + 7] > b[9*8 + 7]) << 7;
            
            f[10] |= (unsigned char) (a[10*8 + 0] > b[10*8 + 0]) << 0;
            f[10] |= (unsigned char) (a[10*8 + 1] > b[10*8 + 1]) << 1;
            f[10] |= (unsigned char) (a[10*8 + 2] > b[10*8 + 2]) << 2;
            f[10] |= (unsigned char) (a[10*8 + 3] > b[10*8 + 3]) << 3;
            f[10] |= (unsigned char) (a[10*8 + 4] > b[10*8 + 4]) << 4;
            f[10] |= (unsigned char) (a[10*8 + 5] > b[10*8 + 5]) << 5;
            f[10] |= (unsigned char) (a[10*8 + 6] > b[10*8 + 6]) << 6;
            f[10] |= (unsigned char) (a[10*8 + 7] > b[10*8 + 7]) << 7;
            
            f[11] |= (unsigned char) (a[11*8 + 0] > b[11*8 + 0]) << 0;
            f[11] |= (unsigned char) (a[11*8 + 1] > b[11*8 + 1]) << 1;
            f[11] |= (unsigned char) (a[11*8 + 2] > b[11*8 + 2]) << 2;
            f[11] |= (unsigned char) (a[11*8 + 3] > b[11*8 + 3]) << 3;
            f[11] |= (unsigned char) (a[11*8 + 4] > b[11*8 + 4]) << 4;
            f[11] |= (unsigned char) (a[11*8 + 5] > b[11*8 + 5]) << 5;
            f[11] |= (unsigned char) (a[11*8 + 6] > b[11*8 + 6]) << 6;
            f[11] |= (unsigned char) (a[11*8 + 7] > b[11*8 + 7]) << 7;
            
            f[12] |= (unsigned char) (a[12*8 + 0] > b[12*8 + 0]) << 0;
            f[12] |= (unsigned char) (a[12*8 + 1] > b[12*8 + 1]) << 1;
            f[12] |= (unsigned char) (a[12*8 + 2] > b[12*8 + 2]) << 2;
            f[12] |= (unsigned char) (a[12*8 + 3] > b[12*8 + 3]) << 3;
            f[12] |= (unsigned char) (a[12*8 + 4] > b[12*8 + 4]) << 4;
            f[12] |= (unsigned char) (a[12*8 + 5] > b[12*8 + 5]) << 5;
            f[12] |= (unsigned char) (a[12*8 + 6] > b[12*8 + 6]) << 6;
            f[12] |= (unsigned char) (a[12*8 + 7] > b[12*8 + 7]) << 7;
            
            f[13] |= (unsigned char) (a[13*8 + 0] > b[13*8 + 0]) << 0;
            f[13] |= (unsigned char) (a[13*8 + 1] > b[13*8 + 1]) << 1;
            f[13] |= (unsigned char) (a[13*8 + 2] > b[13*8 + 2]) << 2;
            f[13] |= (unsigned char) (a[13*8 + 3] > b[13*8 + 3]) << 3;
            f[13] |= (unsigned char) (a[13*8 + 4] > b[13*8 + 4]) << 4;
            f[13] |= (unsigned char) (a[13*8 + 5] > b[13*8 + 5]) << 5;
            f[13] |= (unsigned char) (a[13*8 + 6] > b[13*8 + 6]) << 6;
            f[13] |= (unsigned char) (a[13*8 + 7] > b[13*8 + 7]) << 7;
            
            f[14] |= (unsigned char) (a[14*8 + 0] > b[14*8 + 0]) << 0;
            f[14] |= (unsigned char) (a[14*8 + 1] > b[14*8 + 1]) << 1;
            f[14] |= (unsigned char) (a[14*8 + 2] > b[14*8 + 2]) << 2;
            f[14] |= (unsigned char) (a[14*8 + 3] > b[14*8 + 3]) << 3;
            f[14] |= (unsigned char) (a[14*8 + 4] > b[14*8 + 4]) << 4;
            f[14] |= (unsigned char) (a[14*8 + 5] > b[14*8 + 5]) << 5;
            f[14] |= (unsigned char) (a[14*8 + 6] > b[14*8 + 6]) << 6;
            f[14] |= (unsigned char) (a[14*8 + 7] > b[14*8 + 7]) << 7;
            
            f[15] |= (unsigned char) (a[15*8 + 0] > b[15*8 + 0]) << 0;
            f[15] |= (unsigned char) (a[15*8 + 1] > b[15*8 + 1]) << 1;
            f[15] |= (unsigned char) (a[15*8 + 2] > b[15*8 + 2]) << 2;
            f[15] |= (unsigned char) (a[15*8 + 3] > b[15*8 + 3]) << 3;
            f[15] |= (unsigned char) (a[15*8 + 4] > b[15*8 + 4]) << 4;
            f[15] |= (unsigned char) (a[15*8 + 5] > b[15*8 + 5]) << 5;
            f[15] |= (unsigned char) (a[15*8 + 6] > b[15*8 + 6]) << 6;
            f[15] |= (unsigned char) (a[15*8 + 7] > b[15*8 + 7]) << 7;
            
            f[16] |= (unsigned char) (a[16*8 + 0] > b[16*8 + 0]) << 0;
            f[16] |= (unsigned char) (a[16*8 + 1] > b[16*8 + 1]) << 1;
            f[16] |= (unsigned char) (a[16*8 + 2] > b[16*8 + 2]) << 2;
            f[16] |= (unsigned char) (a[16*8 + 3] > b[16*8 + 3]) << 3;
            f[16] |= (unsigned char) (a[16*8 + 4] > b[16*8 + 4]) << 4;
            f[16] |= (unsigned char) (a[16*8 + 5] > b[16*8 + 5]) << 5;
            f[16] |= (unsigned char) (a[16*8 + 6] > b[16*8 + 6]) << 6;
            f[16] |= (unsigned char) (a[16*8 + 7] > b[16*8 + 7]) << 7;
            
            f[17] |= (unsigned char) (a[17*8 + 0] > b[17*8 + 0]) << 0;
            f[17] |= (unsigned char) (a[17*8 + 1] > b[17*8 + 1]) << 1;
            f[17] |= (unsigned char) (a[17*8 + 2] > b[17*8 + 2]) << 2;
            f[17] |= (unsigned char) (a[17*8 + 3] > b[17*8 + 3]) << 3;
            f[17] |= (unsigned char) (a[17*8 + 4] > b[17*8 + 4]) << 4;
            f[17] |= (unsigned char) (a[17*8 + 5] > b[17*8 + 5]) << 5;
            f[17] |= (unsigned char) (a[17*8 + 6] > b[17*8 + 6]) << 6;
            f[17] |= (unsigned char) (a[17*8 + 7] > b[17*8 + 7]) << 7;
            
            f[18] |= (unsigned char) (a[18*8 + 0] > b[18*8 + 0]) << 0;
            f[18] |= (unsigned char) (a[18*8 + 1] > b[18*8 + 1]) << 1;
            f[18] |= (unsigned char) (a[18*8 + 2] > b[18*8 + 2]) << 2;
            f[18] |= (unsigned char) (a[18*8 + 3] > b[18*8 + 3]) << 3;
            f[18] |= (unsigned char) (a[18*8 + 4] > b[18*8 + 4]) << 4;
            f[18] |= (unsigned char) (a[18*8 + 5] > b[18*8 + 5]) << 5;
            f[18] |= (unsigned char) (a[18*8 + 6] > b[18*8 + 6]) << 6;
            f[18] |= (unsigned char) (a[18*8 + 7] > b[18*8 + 7]) << 7;
            
            f[19] |= (unsigned char) (a[19*8 + 0] > b[19*8 + 0]) << 0;
            f[19] |= (unsigned char) (a[19*8 + 1] > b[19*8 + 1]) << 1;
            f[19] |= (unsigned char) (a[19*8 + 2] > b[19*8 + 2]) << 2;
            f[19] |= (unsigned char) (a[19*8 + 3] > b[19*8 + 3]) << 3;
            f[19] |= (unsigned char) (a[19*8 + 4] > b[19*8 + 4]) << 4;
            f[19] |= (unsigned char) (a[19*8 + 5] > b[19*8 + 5]) << 5;
            f[19] |= (unsigned char) (a[19*8 + 6] > b[19*8 + 6]) << 6;
            f[19] |= (unsigned char) (a[19*8 + 7] > b[19*8 + 7]) << 7;
            
            f[20] |= (unsigned char) (a[20*8 + 0] > b[20*8 + 0]) << 0;
            f[20] |= (unsigned char) (a[20*8 + 1] > b[20*8 + 1]) << 1;
            f[20] |= (unsigned char) (a[20*8 + 2] > b[20*8 + 2]) << 2;
            f[20] |= (unsigned char) (a[20*8 + 3] > b[20*8 + 3]) << 3;
            f[20] |= (unsigned char) (a[20*8 + 4] > b[20*8 + 4]) << 4;
            f[20] |= (unsigned char) (a[20*8 + 5] > b[20*8 + 5]) << 5;
            f[20] |= (unsigned char) (a[20*8 + 6] > b[20*8 + 6]) << 6;
            f[20] |= (unsigned char) (a[20*8 + 7] > b[20*8 + 7]) << 7;
            
            f[21] |= (unsigned char) (a[21*8 + 0] > b[21*8 + 0]) << 0;
            f[21] |= (unsigned char) (a[21*8 + 1] > b[21*8 + 1]) << 1;
            f[21] |= (unsigned char) (a[21*8 + 2] > b[21*8 + 2]) << 2;
            f[21] |= (unsigned char) (a[21*8 + 3] > b[21*8 + 3]) << 3;
            f[21] |= (unsigned char) (a[21*8 + 4] > b[21*8 + 4]) << 4;
            f[21] |= (unsigned char) (a[21*8 + 5] > b[21*8 + 5]) << 5;
            f[21] |= (unsigned char) (a[21*8 + 6] > b[21*8 + 6]) << 6;
            f[21] |= (unsigned char) (a[21*8 + 7] > b[21*8 + 7]) << 7;
            
            f[22] |= (unsigned char) (a[22*8 + 0] > b[22*8 + 0]) << 0;
            f[22] |= (unsigned char) (a[22*8 + 1] > b[22*8 + 1]) << 1;
            f[22] |= (unsigned char) (a[22*8 + 2] > b[22*8 + 2]) << 2;
            f[22] |= (unsigned char) (a[22*8 + 3] > b[22*8 + 3]) << 3;
            f[22] |= (unsigned char) (a[22*8 + 4] > b[22*8 + 4]) << 4;
            f[22] |= (unsigned char) (a[22*8 + 5] > b[22*8 + 5]) << 5;
            f[22] |= (unsigned char) (a[22*8 + 6] > b[22*8 + 6]) << 6;
            f[22] |= (unsigned char) (a[22*8 + 7] > b[22*8 + 7]) << 7;
            
            f[23] |= (unsigned char) (a[23*8 + 0] > b[23*8 + 0]) << 0;
            f[23] |= (unsigned char) (a[23*8 + 1] > b[23*8 + 1]) << 1;
            f[23] |= (unsigned char) (a[23*8 + 2] > b[23*8 + 2]) << 2;
            f[23] |= (unsigned char) (a[23*8 + 3] > b[23*8 + 3]) << 3;
            f[23] |= (unsigned char) (a[23*8 + 4] > b[23*8 + 4]) << 4;
            f[23] |= (unsigned char) (a[23*8 + 5] > b[23*8 + 5]) << 5;
            f[23] |= (unsigned char) (a[23*8 + 6] > b[23*8 + 6]) << 6;
            f[23] |= (unsigned char) (a[23*8 + 7] > b[23*8 + 7]) << 7;
            
            f[24] |= (unsigned char) (a[24*8 + 0] > b[24*8 + 0]) << 0;
            f[24] |= (unsigned char) (a[24*8 + 1] > b[24*8 + 1]) << 1;
            f[24] |= (unsigned char) (a[24*8 + 2] > b[24*8 + 2]) << 2;
            f[24] |= (unsigned char) (a[24*8 + 3] > b[24*8 + 3]) << 3;
            f[24] |= (unsigned char) (a[24*8 + 4] > b[24*8 + 4]) << 4;
            f[24] |= (unsigned char) (a[24*8 + 5] > b[24*8 + 5]) << 5;
            f[24] |= (unsigned char) (a[24*8 + 6] > b[24*8 + 6]) << 6;
            f[24] |= (unsigned char) (a[24*8 + 7] > b[24*8 + 7]) << 7;
            
            f[25] |= (unsigned char) (a[25*8 + 0] > b[25*8 + 0]) << 0;
            f[25] |= (unsigned char) (a[25*8 + 1] > b[25*8 + 1]) << 1;
            f[25] |= (unsigned char) (a[25*8 + 2] > b[25*8 + 2]) << 2;
            f[25] |= (unsigned char) (a[25*8 + 3] > b[25*8 + 3]) << 3;
            f[25] |= (unsigned char) (a[25*8 + 4] > b[25*8 + 4]) << 4;
            f[25] |= (unsigned char) (a[25*8 + 5] > b[25*8 + 5]) << 5;
            f[25] |= (unsigned char) (a[25*8 + 6] > b[25*8 + 6]) << 6;
            f[25] |= (unsigned char) (a[25*8 + 7] > b[25*8 + 7]) << 7;
            
            f[26] |= (unsigned char) (a[26*8 + 0] > b[26*8 + 0]) << 0;
            f[26] |= (unsigned char) (a[26*8 + 1] > b[26*8 + 1]) << 1;
            f[26] |= (unsigned char) (a[26*8 + 2] > b[26*8 + 2]) << 2;
            f[26] |= (unsigned char) (a[26*8 + 3] > b[26*8 + 3]) << 3;
            f[26] |= (unsigned char) (a[26*8 + 4] > b[26*8 + 4]) << 4;
            f[26] |= (unsigned char) (a[26*8 + 5] > b[26*8 + 5]) << 5;
            f[26] |= (unsigned char) (a[26*8 + 6] > b[26*8 + 6]) << 6;
            f[26] |= (unsigned char) (a[26*8 + 7] > b[26*8 + 7]) << 7;
            
            f[27] |= (unsigned char) (a[27*8 + 0] > b[27*8 + 0]) << 0;
            f[27] |= (unsigned char) (a[27*8 + 1] > b[27*8 + 1]) << 1;
            f[27] |= (unsigned char) (a[27*8 + 2] > b[27*8 + 2]) << 2;
            f[27] |= (unsigned char) (a[27*8 + 3] > b[27*8 + 3]) << 3;
            f[27] |= (unsigned char) (a[27*8 + 4] > b[27*8 + 4]) << 4;
            f[27] |= (unsigned char) (a[27*8 + 5] > b[27*8 + 5]) << 5;
            f[27] |= (unsigned char) (a[27*8 + 6] > b[27*8 + 6]) << 6;
            f[27] |= (unsigned char) (a[27*8 + 7] > b[27*8 + 7]) << 7;
            
            f[28] |= (unsigned char) (a[28*8 + 0] > b[28*8 + 0]) << 0;
            f[28] |= (unsigned char) (a[28*8 + 1] > b[28*8 + 1]) << 1;
            f[28] |= (unsigned char) (a[28*8 + 2] > b[28*8 + 2]) << 2;
            f[28] |= (unsigned char) (a[28*8 + 3] > b[28*8 + 3]) << 3;
            f[28] |= (unsigned char) (a[28*8 + 4] > b[28*8 + 4]) << 4;
            f[28] |= (unsigned char) (a[28*8 + 5] > b[28*8 + 5]) << 5;
            f[28] |= (unsigned char) (a[28*8 + 6] > b[28*8 + 6]) << 6;
            f[28] |= (unsigned char) (a[28*8 + 7] > b[28*8 + 7]) << 7;
            
            f[29] |= (unsigned char) (a[29*8 + 0] > b[29*8 + 0]) << 0;
            f[29] |= (unsigned char) (a[29*8 + 1] > b[29*8 + 1]) << 1;
            f[29] |= (unsigned char) (a[29*8 + 2] > b[29*8 + 2]) << 2;
            f[29] |= (unsigned char) (a[29*8 + 3] > b[29*8 + 3]) << 3;
            f[29] |= (unsigned char) (a[29*8 + 4] > b[29*8 + 4]) << 4;
            f[29] |= (unsigned char) (a[29*8 + 5] > b[29*8 + 5]) << 5;
            f[29] |= (unsigned char) (a[29*8 + 6] > b[29*8 + 6]) << 6;
            f[29] |= (unsigned char) (a[29*8 + 7] > b[29*8 + 7]) << 7;
            
            f[30] |= (unsigned char) (a[30*8 + 0] > b[30*8 + 0]) << 0;
            f[30] |= (unsigned char) (a[30*8 + 1] > b[30*8 + 1]) << 1;
            f[30] |= (unsigned char) (a[30*8 + 2] > b[30*8 + 2]) << 2;
            f[30] |= (unsigned char) (a[30*8 + 3] > b[30*8 + 3]) << 3;
            f[30] |= (unsigned char) (a[30*8 + 4] > b[30*8 + 4]) << 4;
            f[30] |= (unsigned char) (a[30*8 + 5] > b[30*8 + 5]) << 5;
            f[30] |= (unsigned char) (a[30*8 + 6] > b[30*8 + 6]) << 6;
            f[30] |= (unsigned char) (a[30*8 + 7] > b[30*8 + 7]) << 7;
            
            f[31] |= (unsigned char) (a[31*8 + 0] > b[31*8 + 0]) << 0;
            f[31] |= (unsigned char) (a[31*8 + 1] > b[31*8 + 1]) << 1;
            f[31] |= (unsigned char) (a[31*8 + 2] > b[31*8 + 2]) << 2;
            f[31] |= (unsigned char) (a[31*8 + 3] > b[31*8 + 3]) << 3;
            f[31] |= (unsigned char) (a[31*8 + 4] > b[31*8 + 4]) << 4;
            f[31] |= (unsigned char) (a[31*8 + 5] > b[31*8 + 5]) << 5;
            f[31] |= (unsigned char) (a[31*8 + 6] > b[31*8 + 6]) << 6;
            f[31] |= (unsigned char) (a[31*8 + 7] > b[31*8 + 7]) << 7;
            
            
            _mm_store_si128((__m128i*)(bd + 0*n_rows_bd),_mm_load_si128((const __m128i *)(f + 0*64)));
            _mm_store_si128((__m128i*)(bd + 1*n_rows_bd),_mm_load_si128((const __m128i *)(f + 1*64)));
            _mm_store_si128((__m128i*)(bd + 2*n_rows_bd),_mm_load_si128((const __m128i *)(f + 2*64)));
            _mm_store_si128((__m128i*)(bd + 3*n_rows_bd),_mm_load_si128((const __m128i *)(f + 3*64)));
            

        }
    }
}


int BRIEF::diag_length_pattern = 17;
int BRIEF::gaussian_bit_pattern_31_x_a[256] = { 8,4,-11,7,2,1,-2,-13,-13,10,-13,-11,7,-4,-13,-9,12,-3,-6,11,4,5,3,-8,-2,-13,-7,-4,-10,5,5,1,9,4,2,-4,-8,4,0,-13,-3,-6,8,0,7,-13,10,-6,10,-13,-13,3,5,-1,3,2,-13,-13,-13,-7,6,-9,-2,-12,3,-7,-3,2,-11,-1,5,-4,-9,-12,10,7,-7,-4,7,-7,-13,-3,7,-13,1,2,-4,-1,7,1,9,-1,-13,7,12,6,5,2,3,2,9,-8,-11,1,6,2,6,3,7,-11,-10,-5,-10,8,4,-10,4,-2,-5,7,-9,-5,8,-9,1,7,-2,11,-12,3,5,0,-9,0,-1,5,3,-13,-5,-4,6,-7,-13,1,4,-2,2,-2,4,-6,-3,7,4,-13,7,7,-7,-8,-13,2,10,-6,8,2,-11,-12,-11,5,-2,-1,-13,-10,-3,2,-9,-4,-4,-6,6,-13,11,7,-1,-4,-7,-13,-7,-8,-5,-13,1,1,9,5,-1,-9,-1,-13,8,2,7,-10,-10,4,3,-4,5,4,-9,0,-12,3,-10,8,-8,2,10,6,-7,-3,-1,-3,-8,4,2,6,3,11,-3,4,2,-10,-13,-13,6,0,-13,-9,-13,5,2,-1,9,11,3,-1,3,-13,5,8,7,-10,7,9,7,-1};
int BRIEF::gaussian_bit_pattern_31_y_a[256] = {-3,2,9,-12,-13,-7,-10,-13,-3,4,-8,7,7,-5,2,0,-6,6,-13,-13,7,-3,-7,-7,11,12,3,2,-12,-12,-6,0,11,7,-1,-12,-5,11,-8,-2,-2,9,12,9,-5,-6,7,-3,-9,8,0,3,7,7,-10,-4,0,-7,3,12,-10,-1,-5,5,-10,-7,-2,9,-13,6,-3,-13,-6,-10,2,12,-13,9,-1,6,11,7,-8,-7,-3,-6,3,-13,1,-1,1,-9,-13,7,-5,3,-13,-12,8,6,-12,4,12,12,-9,3,3,-3,8,-5,11,-8,5,-1,-6,12,-2,0,-8,-6,-13,-13,-8,-11,-8,-4,1,-6,-9,7,5,-4,12,7,2,11,5,-4,9,-7,5,6,6,-10,1,-2,-12,-13,1,-10,-13,5,-2,9,1,-8,-4,11,6,4,-5,-5,-3,-12,-2,-13,0,-3,-13,-8,-11,-2,9,-3,-13,6,12,-11,-3,11,11,-5,12,-8,1,-12,-2,5,-1,7,5,0,12,-8,11,-3,-10,1,-11,-13,-13,-10,-8,-6,12,2,-13,-13,9,3,1,2,-10,-13,-12,2,6,8,10,-9,-13,-7,-2,2,-5,-9,-1,-1,0,-11,-4,-6,7,12,0,-1,3,8,-6,-9,7,-6,5,-3,0,4,-6,0,8,9,-4,4,3,-7,0,-6};
int BRIEF::gaussian_bit_pattern_31_x_b[256] = {9,7,-8,12,2,1,-2,-11,-12,11,-8,-9,12,-3,-12,-7,12,-2,-4,12,5,10,6,-6,-1,-8,-5,-3,-6,6,7,4,11,4,4,-2,-7,9,1,-8,-2,-4,10,1,11,-11,12,-6,12,-8,-8,7,10,1,5,3,-13,-12,-11,-4,12,-7,0,-7,8,-4,-1,5,-5,0,5,-4,-9,-8,12,12,-6,-3,12,-5,-12,-2,12,-11,12,3,-2,1,8,3,12,-1,-10,10,12,7,6,2,4,12,10,-7,-4,2,7,3,11,8,9,-6,-5,-3,-9,12,6,-8,6,-2,-5,10,-8,-5,9,-9,1,9,-1,12,-6,7,10,2,-5,2,1,7,6,-8,-3,-3,8,-6,-5,3,8,2,12,0,9,-3,-1,12,5,-9,8,7,-7,-7,-12,3,12,-6,9,2,-10,-7,-10,11,-1,0,-12,-10,-2,3,-4,-3,-2,-4,6,-5,12,12,0,-3,-6,-8,-6,-6,-4,-8,5,10,10,10,1,-6,1,-8,10,3,12,-5,-8,8,8,-3,10,5,-4,3,-6,4,-10,12,-6,3,11,8,-6,-3,-1,-3,-8,12,3,11,7,12,-3,4,2,-8,-11,-11,11,1,-9,-6,-8,8,3,-1,11,12,3,0,4,-10,12,9,8,-10,12,10,12,0};
int BRIEF::gaussian_bit_pattern_31_y_b[256] = {5,-12,2,-13,12,6,-4,-8,-9,9,-9,12,6,0,-3,5,-1,12,-8,-8,1,-3,12,-2,-10,10,-3,7,11,-7,-1,-5,-13,12,4,7,-10,12,-13,2,3,-9,7,3,-10,0,1,12,-4,-12,-4,8,-7,-12,6,-10,5,12,8,7,8,-6,12,5,-13,5,-7,-11,-13,-1,2,12,6,-4,-3,12,5,4,2,1,5,-6,-7,-12,12,0,-13,9,-6,12,6,3,5,12,9,11,10,3,-6,-13,3,9,-6,-8,-4,-2,0,-8,3,-4,10,12,0,-6,-11,7,7,12,2,12,-8,-2,-13,0,-2,1,-4,-11,4,12,8,8,-13,12,7,-9,-8,9,-3,-12,0,12,-2,10,-4,-13,12,-6,3,-5,1,-11,-7,-5,6,6,1,-8,-8,9,3,7,-8,8,3,-9,-5,8,12,9,-5,11,-13,2,0,-10,-7,9,11,5,6,-2,7,-2,7,-13,-8,-9,5,10,-13,-13,-1,-9,-13,2,12,-10,-6,-6,-9,-7,-13,5,-13,-3,-12,-1,3,-9,1,-8,9,12,-5,7,-8,-12,5,9,5,4,3,12,11,-13,12,4,6,12,1,1,1,-13,-13,4,-2,-3,-2,10,-9,-1,-2,-8,5,10,5,5,11,-6,-12,9,4,-2,-2,-11};