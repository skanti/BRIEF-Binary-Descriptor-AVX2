

#include "Types.h"
#include "BRIEF.h"
#include <cmath>
#include <iostream>
#include "immintrin.h"


extern "C" void kernel1(int *p_x_a, int *p_y_a, int *p_x_b, int *p_y_b,
                               int *i_x_a, int *i_y_a, int *i_x_b, int *i_y_b,
                               float cos, float sin);

BRIEF::BRIEF(int n_rows_, int n_cols_)
    : n_rows(n_rows_), n_cols(n_cols_), bd(n_rows, n_cols) {};


void
BRIEF::rbrief(unsigned char *image_src, const int height_image, const int width_image, const int n_channels,
             const int stride_image, const int *x, const int *y, const float *angle, const int n_features) {
    for (int j = 0; j < n_features; j++) {
        if ((x[j] > diag_length_pattern) && x[j] < (width_image - diag_length_pattern)
            && (y[j] > diag_length_pattern) && y[j] < (height_image - diag_length_pattern)) {
            float cos_angle = std::cos(angle[j]);
            float sin_angle = std::sin(angle[j]);
            float tri0[4]  __attribute__((aligned(32))) = {cos_angle, sin_angle, cos_angle, sin_angle};
            float tri1[4]  __attribute__((aligned(32))) = {sin_angle,cos_angle, sin_angle, cos_angle};
            unsigned char *image_center = image_src + y[j] * stride_image + x[j] * n_channels;
            // N_DIM_BINARYDESCRIPTOR / SIZE_BITS_HAMING = 4
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

            int32_t f[8] __attribute__((aligned(32)));
            f[0] |= (*(image_center + i_y_a[0*32 + 0]*stride_image + i_x_a[0*32 + 0])
                        > *(image_center + i_y_b[0*32 + 0]*stride_image + i_x_b[0*32 + 0])) << 0;
            f[0] |= (*(image_center + i_y_a[0*32 + 1]*stride_image + i_x_a[0*32 + 1])
                        > *(image_center + i_y_b[0*32 + 1]*stride_image + i_x_b[0*32 + 1])) << 1;
            f[0] |= (*(image_center + i_y_a[0*32 + 2]*stride_image + i_x_a[0*32 + 2])
                        > *(image_center + i_y_b[0*32 + 2]*stride_image + i_x_b[0*32 + 2])) << 2;
            f[0] |= (*(image_center + i_y_a[0*32 + 3]*stride_image + i_x_a[0*32 + 3])
                        > *(image_center + i_y_b[0*32 + 3]*stride_image + i_x_b[0*32 + 3])) << 3;
            f[0] |= (*(image_center + i_y_a[0*32 + 4]*stride_image + i_x_a[0*32 + 4])
                        > *(image_center + i_y_b[0*32 + 4]*stride_image + i_x_b[0*32 + 4])) << 4;
            f[0] |= (*(image_center + i_y_a[0*32 + 5]*stride_image + i_x_a[0*32 + 5])
                        > *(image_center + i_y_b[0*32 + 5]*stride_image + i_x_b[0*32 + 5])) << 5;
            f[0] |= (*(image_center + i_y_a[0*32 + 6]*stride_image + i_x_a[0*32 + 6])
                        > *(image_center + i_y_b[0*32 + 6]*stride_image + i_x_b[0*32 + 6])) << 6;
            f[0] |= (*(image_center + i_y_a[0*32 + 7]*stride_image + i_x_a[0*32 + 7])
                        > *(image_center + i_y_b[0*32 + 7]*stride_image + i_x_b[0*32 + 7])) << 7;
            f[0] |= (*(image_center + i_y_a[0*32 + 8]*stride_image + i_x_a[0*32 + 8])
                        > *(image_center + i_y_b[0*32 + 8]*stride_image + i_x_b[0*32 + 8])) << 8;
            f[0] |= (*(image_center + i_y_a[0*32 + 9]*stride_image + i_x_a[0*32 + 9])
                        > *(image_center + i_y_b[0*32 + 9]*stride_image + i_x_b[0*32 + 9])) << 9;
            f[0] |= (*(image_center + i_y_a[0*32 + 10]*stride_image + i_x_a[0*32 + 10])
                        > *(image_center + i_y_b[0*32 + 10]*stride_image + i_x_b[0*32 + 10])) << 10;
            f[0] |= (*(image_center + i_y_a[0*32 + 11]*stride_image + i_x_a[0*32 + 11])
                        > *(image_center + i_y_b[0*32 + 11]*stride_image + i_x_b[0*32 + 11])) << 11;
            f[0] |= (*(image_center + i_y_a[0*32 + 12]*stride_image + i_x_a[0*32 + 12])
                        > *(image_center + i_y_b[0*32 + 12]*stride_image + i_x_b[0*32 + 12])) << 12;
            f[0] |= (*(image_center + i_y_a[0*32 + 13]*stride_image + i_x_a[0*32 + 13])
                        > *(image_center + i_y_b[0*32 + 13]*stride_image + i_x_b[0*32 + 13])) << 13;
            f[0] |= (*(image_center + i_y_a[0*32 + 14]*stride_image + i_x_a[0*32 + 14])
                        > *(image_center + i_y_b[0*32 + 14]*stride_image + i_x_b[0*32 + 14])) << 14;
            f[0] |= (*(image_center + i_y_a[0*32 + 15]*stride_image + i_x_a[0*32 + 15])
                        > *(image_center + i_y_b[0*32 + 15]*stride_image + i_x_b[0*32 + 15])) << 15;
            f[0] |= (*(image_center + i_y_a[0*32 + 16]*stride_image + i_x_a[0*32 + 16])
                        > *(image_center + i_y_b[0*32 + 16]*stride_image + i_x_b[0*32 + 16])) << 16;
            f[0] |= (*(image_center + i_y_a[0*32 + 17]*stride_image + i_x_a[0*32 + 17])
                        > *(image_center + i_y_b[0*32 + 17]*stride_image + i_x_b[0*32 + 17])) << 17;
            f[0] |= (*(image_center + i_y_a[0*32 + 18]*stride_image + i_x_a[0*32 + 18])
                        > *(image_center + i_y_b[0*32 + 18]*stride_image + i_x_b[0*32 + 18])) << 18;
            f[0] |= (*(image_center + i_y_a[0*32 + 19]*stride_image + i_x_a[0*32 + 19])
                        > *(image_center + i_y_b[0*32 + 19]*stride_image + i_x_b[0*32 + 19])) << 19;
            f[0] |= (*(image_center + i_y_a[0*32 + 20]*stride_image + i_x_a[0*32 + 20])
                        > *(image_center + i_y_b[0*32 + 20]*stride_image + i_x_b[0*32 + 20])) << 20;
            f[0] |= (*(image_center + i_y_a[0*32 + 21]*stride_image + i_x_a[0*32 + 21])
                        > *(image_center + i_y_b[0*32 + 21]*stride_image + i_x_b[0*32 + 21])) << 21;
            f[0] |= (*(image_center + i_y_a[0*32 + 22]*stride_image + i_x_a[0*32 + 22])
                        > *(image_center + i_y_b[0*32 + 22]*stride_image + i_x_b[0*32 + 22])) << 22;
            f[0] |= (*(image_center + i_y_a[0*32 + 23]*stride_image + i_x_a[0*32 + 23])
                        > *(image_center + i_y_b[0*32 + 23]*stride_image + i_x_b[0*32 + 23])) << 23;
            f[0] |= (*(image_center + i_y_a[0*32 + 24]*stride_image + i_x_a[0*32 + 24])
                        > *(image_center + i_y_b[0*32 + 24]*stride_image + i_x_b[0*32 + 24])) << 24;
            f[0] |= (*(image_center + i_y_a[0*32 + 25]*stride_image + i_x_a[0*32 + 25])
                        > *(image_center + i_y_b[0*32 + 25]*stride_image + i_x_b[0*32 + 25])) << 25;
            f[0] |= (*(image_center + i_y_a[0*32 + 26]*stride_image + i_x_a[0*32 + 26])
                        > *(image_center + i_y_b[0*32 + 26]*stride_image + i_x_b[0*32 + 26])) << 26;
            f[0] |= (*(image_center + i_y_a[0*32 + 27]*stride_image + i_x_a[0*32 + 27])
                        > *(image_center + i_y_b[0*32 + 27]*stride_image + i_x_b[0*32 + 27])) << 27;
            f[0] |= (*(image_center + i_y_a[0*32 + 28]*stride_image + i_x_a[0*32 + 28])
                        > *(image_center + i_y_b[0*32 + 28]*stride_image + i_x_b[0*32 + 28])) << 28;
            f[0] |= (*(image_center + i_y_a[0*32 + 29]*stride_image + i_x_a[0*32 + 29])
                        > *(image_center + i_y_b[0*32 + 29]*stride_image + i_x_b[0*32 + 29])) << 29;
            f[0] |= (*(image_center + i_y_a[0*32 + 30]*stride_image + i_x_a[0*32 + 30])
                        > *(image_center + i_y_b[0*32 + 30]*stride_image + i_x_b[0*32 + 30])) << 30;
            f[0] |= (*(image_center + i_y_a[0*32 + 31]*stride_image + i_x_a[0*32 + 31])
                        > *(image_center + i_y_b[0*32 + 31]*stride_image + i_x_b[0*32 + 31])) << 31;
            
            f[1] |= (*(image_center + i_y_a[1*32 + 0]*stride_image + i_x_a[1*32 + 0])
                        > *(image_center + i_y_b[1*32 + 0]*stride_image + i_x_b[1*32 + 0])) << 0;
            f[1] |= (*(image_center + i_y_a[1*32 + 1]*stride_image + i_x_a[1*32 + 1])
                        > *(image_center + i_y_b[1*32 + 1]*stride_image + i_x_b[1*32 + 1])) << 1;
            f[1] |= (*(image_center + i_y_a[1*32 + 2]*stride_image + i_x_a[1*32 + 2])
                        > *(image_center + i_y_b[1*32 + 2]*stride_image + i_x_b[1*32 + 2])) << 2;
            f[1] |= (*(image_center + i_y_a[1*32 + 3]*stride_image + i_x_a[1*32 + 3])
                        > *(image_center + i_y_b[1*32 + 3]*stride_image + i_x_b[1*32 + 3])) << 3;
            f[1] |= (*(image_center + i_y_a[1*32 + 4]*stride_image + i_x_a[1*32 + 4])
                        > *(image_center + i_y_b[1*32 + 4]*stride_image + i_x_b[1*32 + 4])) << 4;
            f[1] |= (*(image_center + i_y_a[1*32 + 5]*stride_image + i_x_a[1*32 + 5])
                        > *(image_center + i_y_b[1*32 + 5]*stride_image + i_x_b[1*32 + 5])) << 5;
            f[1] |= (*(image_center + i_y_a[1*32 + 6]*stride_image + i_x_a[1*32 + 6])
                        > *(image_center + i_y_b[1*32 + 6]*stride_image + i_x_b[1*32 + 6])) << 6;
            f[1] |= (*(image_center + i_y_a[1*32 + 7]*stride_image + i_x_a[1*32 + 7])
                        > *(image_center + i_y_b[1*32 + 7]*stride_image + i_x_b[1*32 + 7])) << 7;
            f[1] |= (*(image_center + i_y_a[1*32 + 8]*stride_image + i_x_a[1*32 + 8])
                        > *(image_center + i_y_b[1*32 + 8]*stride_image + i_x_b[1*32 + 8])) << 8;
            f[1] |= (*(image_center + i_y_a[1*32 + 9]*stride_image + i_x_a[1*32 + 9])
                        > *(image_center + i_y_b[1*32 + 9]*stride_image + i_x_b[1*32 + 9])) << 9;
            f[1] |= (*(image_center + i_y_a[1*32 + 10]*stride_image + i_x_a[1*32 + 10])
                        > *(image_center + i_y_b[1*32 + 10]*stride_image + i_x_b[1*32 + 10])) << 10;
            f[1] |= (*(image_center + i_y_a[1*32 + 11]*stride_image + i_x_a[1*32 + 11])
                        > *(image_center + i_y_b[1*32 + 11]*stride_image + i_x_b[1*32 + 11])) << 11;
            f[1] |= (*(image_center + i_y_a[1*32 + 12]*stride_image + i_x_a[1*32 + 12])
                        > *(image_center + i_y_b[1*32 + 12]*stride_image + i_x_b[1*32 + 12])) << 12;
            f[1] |= (*(image_center + i_y_a[1*32 + 13]*stride_image + i_x_a[1*32 + 13])
                        > *(image_center + i_y_b[1*32 + 13]*stride_image + i_x_b[1*32 + 13])) << 13;
            f[1] |= (*(image_center + i_y_a[1*32 + 14]*stride_image + i_x_a[1*32 + 14])
                        > *(image_center + i_y_b[1*32 + 14]*stride_image + i_x_b[1*32 + 14])) << 14;
            f[1] |= (*(image_center + i_y_a[1*32 + 15]*stride_image + i_x_a[1*32 + 15])
                        > *(image_center + i_y_b[1*32 + 15]*stride_image + i_x_b[1*32 + 15])) << 15;
            f[1] |= (*(image_center + i_y_a[1*32 + 16]*stride_image + i_x_a[1*32 + 16])
                        > *(image_center + i_y_b[1*32 + 16]*stride_image + i_x_b[1*32 + 16])) << 16;
            f[1] |= (*(image_center + i_y_a[1*32 + 17]*stride_image + i_x_a[1*32 + 17])
                        > *(image_center + i_y_b[1*32 + 17]*stride_image + i_x_b[1*32 + 17])) << 17;
            f[1] |= (*(image_center + i_y_a[1*32 + 18]*stride_image + i_x_a[1*32 + 18])
                        > *(image_center + i_y_b[1*32 + 18]*stride_image + i_x_b[1*32 + 18])) << 18;
            f[1] |= (*(image_center + i_y_a[1*32 + 19]*stride_image + i_x_a[1*32 + 19])
                        > *(image_center + i_y_b[1*32 + 19]*stride_image + i_x_b[1*32 + 19])) << 19;
            f[1] |= (*(image_center + i_y_a[1*32 + 20]*stride_image + i_x_a[1*32 + 20])
                        > *(image_center + i_y_b[1*32 + 20]*stride_image + i_x_b[1*32 + 20])) << 20;
            f[1] |= (*(image_center + i_y_a[1*32 + 21]*stride_image + i_x_a[1*32 + 21])
                        > *(image_center + i_y_b[1*32 + 21]*stride_image + i_x_b[1*32 + 21])) << 21;
            f[1] |= (*(image_center + i_y_a[1*32 + 22]*stride_image + i_x_a[1*32 + 22])
                        > *(image_center + i_y_b[1*32 + 22]*stride_image + i_x_b[1*32 + 22])) << 22;
            f[1] |= (*(image_center + i_y_a[1*32 + 23]*stride_image + i_x_a[1*32 + 23])
                        > *(image_center + i_y_b[1*32 + 23]*stride_image + i_x_b[1*32 + 23])) << 23;
            f[1] |= (*(image_center + i_y_a[1*32 + 24]*stride_image + i_x_a[1*32 + 24])
                        > *(image_center + i_y_b[1*32 + 24]*stride_image + i_x_b[1*32 + 24])) << 24;
            f[1] |= (*(image_center + i_y_a[1*32 + 25]*stride_image + i_x_a[1*32 + 25])
                        > *(image_center + i_y_b[1*32 + 25]*stride_image + i_x_b[1*32 + 25])) << 25;
            f[1] |= (*(image_center + i_y_a[1*32 + 26]*stride_image + i_x_a[1*32 + 26])
                        > *(image_center + i_y_b[1*32 + 26]*stride_image + i_x_b[1*32 + 26])) << 26;
            f[1] |= (*(image_center + i_y_a[1*32 + 27]*stride_image + i_x_a[1*32 + 27])
                        > *(image_center + i_y_b[1*32 + 27]*stride_image + i_x_b[1*32 + 27])) << 27;
            f[1] |= (*(image_center + i_y_a[1*32 + 28]*stride_image + i_x_a[1*32 + 28])
                        > *(image_center + i_y_b[1*32 + 28]*stride_image + i_x_b[1*32 + 28])) << 28;
            f[1] |= (*(image_center + i_y_a[1*32 + 29]*stride_image + i_x_a[1*32 + 29])
                        > *(image_center + i_y_b[1*32 + 29]*stride_image + i_x_b[1*32 + 29])) << 29;
            f[1] |= (*(image_center + i_y_a[1*32 + 30]*stride_image + i_x_a[1*32 + 30])
                        > *(image_center + i_y_b[1*32 + 30]*stride_image + i_x_b[1*32 + 30])) << 30;
            f[1] |= (*(image_center + i_y_a[1*32 + 31]*stride_image + i_x_a[1*32 + 31])
                        > *(image_center + i_y_b[1*32 + 31]*stride_image + i_x_b[1*32 + 31])) << 31;
            
            f[2] |= (*(image_center + i_y_a[2*32 + 0]*stride_image + i_x_a[2*32 + 0])
                        > *(image_center + i_y_b[2*32 + 0]*stride_image + i_x_b[2*32 + 0])) << 0;
            f[2] |= (*(image_center + i_y_a[2*32 + 1]*stride_image + i_x_a[2*32 + 1])
                        > *(image_center + i_y_b[2*32 + 1]*stride_image + i_x_b[2*32 + 1])) << 1;
            f[2] |= (*(image_center + i_y_a[2*32 + 2]*stride_image + i_x_a[2*32 + 2])
                        > *(image_center + i_y_b[2*32 + 2]*stride_image + i_x_b[2*32 + 2])) << 2;
            f[2] |= (*(image_center + i_y_a[2*32 + 3]*stride_image + i_x_a[2*32 + 3])
                        > *(image_center + i_y_b[2*32 + 3]*stride_image + i_x_b[2*32 + 3])) << 3;
            f[2] |= (*(image_center + i_y_a[2*32 + 4]*stride_image + i_x_a[2*32 + 4])
                        > *(image_center + i_y_b[2*32 + 4]*stride_image + i_x_b[2*32 + 4])) << 4;
            f[2] |= (*(image_center + i_y_a[2*32 + 5]*stride_image + i_x_a[2*32 + 5])
                        > *(image_center + i_y_b[2*32 + 5]*stride_image + i_x_b[2*32 + 5])) << 5;
            f[2] |= (*(image_center + i_y_a[2*32 + 6]*stride_image + i_x_a[2*32 + 6])
                        > *(image_center + i_y_b[2*32 + 6]*stride_image + i_x_b[2*32 + 6])) << 6;
            f[2] |= (*(image_center + i_y_a[2*32 + 7]*stride_image + i_x_a[2*32 + 7])
                        > *(image_center + i_y_b[2*32 + 7]*stride_image + i_x_b[2*32 + 7])) << 7;
            f[2] |= (*(image_center + i_y_a[2*32 + 8]*stride_image + i_x_a[2*32 + 8])
                        > *(image_center + i_y_b[2*32 + 8]*stride_image + i_x_b[2*32 + 8])) << 8;
            f[2] |= (*(image_center + i_y_a[2*32 + 9]*stride_image + i_x_a[2*32 + 9])
                        > *(image_center + i_y_b[2*32 + 9]*stride_image + i_x_b[2*32 + 9])) << 9;
            f[2] |= (*(image_center + i_y_a[2*32 + 10]*stride_image + i_x_a[2*32 + 10])
                        > *(image_center + i_y_b[2*32 + 10]*stride_image + i_x_b[2*32 + 10])) << 10;
            f[2] |= (*(image_center + i_y_a[2*32 + 11]*stride_image + i_x_a[2*32 + 11])
                        > *(image_center + i_y_b[2*32 + 11]*stride_image + i_x_b[2*32 + 11])) << 11;
            f[2] |= (*(image_center + i_y_a[2*32 + 12]*stride_image + i_x_a[2*32 + 12])
                        > *(image_center + i_y_b[2*32 + 12]*stride_image + i_x_b[2*32 + 12])) << 12;
            f[2] |= (*(image_center + i_y_a[2*32 + 13]*stride_image + i_x_a[2*32 + 13])
                        > *(image_center + i_y_b[2*32 + 13]*stride_image + i_x_b[2*32 + 13])) << 13;
            f[2] |= (*(image_center + i_y_a[2*32 + 14]*stride_image + i_x_a[2*32 + 14])
                        > *(image_center + i_y_b[2*32 + 14]*stride_image + i_x_b[2*32 + 14])) << 14;
            f[2] |= (*(image_center + i_y_a[2*32 + 15]*stride_image + i_x_a[2*32 + 15])
                        > *(image_center + i_y_b[2*32 + 15]*stride_image + i_x_b[2*32 + 15])) << 15;
            f[2] |= (*(image_center + i_y_a[2*32 + 16]*stride_image + i_x_a[2*32 + 16])
                        > *(image_center + i_y_b[2*32 + 16]*stride_image + i_x_b[2*32 + 16])) << 16;
            f[2] |= (*(image_center + i_y_a[2*32 + 17]*stride_image + i_x_a[2*32 + 17])
                        > *(image_center + i_y_b[2*32 + 17]*stride_image + i_x_b[2*32 + 17])) << 17;
            f[2] |= (*(image_center + i_y_a[2*32 + 18]*stride_image + i_x_a[2*32 + 18])
                        > *(image_center + i_y_b[2*32 + 18]*stride_image + i_x_b[2*32 + 18])) << 18;
            f[2] |= (*(image_center + i_y_a[2*32 + 19]*stride_image + i_x_a[2*32 + 19])
                        > *(image_center + i_y_b[2*32 + 19]*stride_image + i_x_b[2*32 + 19])) << 19;
            f[2] |= (*(image_center + i_y_a[2*32 + 20]*stride_image + i_x_a[2*32 + 20])
                        > *(image_center + i_y_b[2*32 + 20]*stride_image + i_x_b[2*32 + 20])) << 20;
            f[2] |= (*(image_center + i_y_a[2*32 + 21]*stride_image + i_x_a[2*32 + 21])
                        > *(image_center + i_y_b[2*32 + 21]*stride_image + i_x_b[2*32 + 21])) << 21;
            f[2] |= (*(image_center + i_y_a[2*32 + 22]*stride_image + i_x_a[2*32 + 22])
                        > *(image_center + i_y_b[2*32 + 22]*stride_image + i_x_b[2*32 + 22])) << 22;
            f[2] |= (*(image_center + i_y_a[2*32 + 23]*stride_image + i_x_a[2*32 + 23])
                        > *(image_center + i_y_b[2*32 + 23]*stride_image + i_x_b[2*32 + 23])) << 23;
            f[2] |= (*(image_center + i_y_a[2*32 + 24]*stride_image + i_x_a[2*32 + 24])
                        > *(image_center + i_y_b[2*32 + 24]*stride_image + i_x_b[2*32 + 24])) << 24;
            f[2] |= (*(image_center + i_y_a[2*32 + 25]*stride_image + i_x_a[2*32 + 25])
                        > *(image_center + i_y_b[2*32 + 25]*stride_image + i_x_b[2*32 + 25])) << 25;
            f[2] |= (*(image_center + i_y_a[2*32 + 26]*stride_image + i_x_a[2*32 + 26])
                        > *(image_center + i_y_b[2*32 + 26]*stride_image + i_x_b[2*32 + 26])) << 26;
            f[2] |= (*(image_center + i_y_a[2*32 + 27]*stride_image + i_x_a[2*32 + 27])
                        > *(image_center + i_y_b[2*32 + 27]*stride_image + i_x_b[2*32 + 27])) << 27;
            f[2] |= (*(image_center + i_y_a[2*32 + 28]*stride_image + i_x_a[2*32 + 28])
                        > *(image_center + i_y_b[2*32 + 28]*stride_image + i_x_b[2*32 + 28])) << 28;
            f[2] |= (*(image_center + i_y_a[2*32 + 29]*stride_image + i_x_a[2*32 + 29])
                        > *(image_center + i_y_b[2*32 + 29]*stride_image + i_x_b[2*32 + 29])) << 29;
            f[2] |= (*(image_center + i_y_a[2*32 + 30]*stride_image + i_x_a[2*32 + 30])
                        > *(image_center + i_y_b[2*32 + 30]*stride_image + i_x_b[2*32 + 30])) << 30;
            f[2] |= (*(image_center + i_y_a[2*32 + 31]*stride_image + i_x_a[2*32 + 31])
                        > *(image_center + i_y_b[2*32 + 31]*stride_image + i_x_b[2*32 + 31])) << 31;
            
            f[3] |= (*(image_center + i_y_a[3*32 + 0]*stride_image + i_x_a[3*32 + 0])
                        > *(image_center + i_y_b[3*32 + 0]*stride_image + i_x_b[3*32 + 0])) << 0;
            f[3] |= (*(image_center + i_y_a[3*32 + 1]*stride_image + i_x_a[3*32 + 1])
                        > *(image_center + i_y_b[3*32 + 1]*stride_image + i_x_b[3*32 + 1])) << 1;
            f[3] |= (*(image_center + i_y_a[3*32 + 2]*stride_image + i_x_a[3*32 + 2])
                        > *(image_center + i_y_b[3*32 + 2]*stride_image + i_x_b[3*32 + 2])) << 2;
            f[3] |= (*(image_center + i_y_a[3*32 + 3]*stride_image + i_x_a[3*32 + 3])
                        > *(image_center + i_y_b[3*32 + 3]*stride_image + i_x_b[3*32 + 3])) << 3;
            f[3] |= (*(image_center + i_y_a[3*32 + 4]*stride_image + i_x_a[3*32 + 4])
                        > *(image_center + i_y_b[3*32 + 4]*stride_image + i_x_b[3*32 + 4])) << 4;
            f[3] |= (*(image_center + i_y_a[3*32 + 5]*stride_image + i_x_a[3*32 + 5])
                        > *(image_center + i_y_b[3*32 + 5]*stride_image + i_x_b[3*32 + 5])) << 5;
            f[3] |= (*(image_center + i_y_a[3*32 + 6]*stride_image + i_x_a[3*32 + 6])
                        > *(image_center + i_y_b[3*32 + 6]*stride_image + i_x_b[3*32 + 6])) << 6;
            f[3] |= (*(image_center + i_y_a[3*32 + 7]*stride_image + i_x_a[3*32 + 7])
                        > *(image_center + i_y_b[3*32 + 7]*stride_image + i_x_b[3*32 + 7])) << 7;
            f[3] |= (*(image_center + i_y_a[3*32 + 8]*stride_image + i_x_a[3*32 + 8])
                        > *(image_center + i_y_b[3*32 + 8]*stride_image + i_x_b[3*32 + 8])) << 8;
            f[3] |= (*(image_center + i_y_a[3*32 + 9]*stride_image + i_x_a[3*32 + 9])
                        > *(image_center + i_y_b[3*32 + 9]*stride_image + i_x_b[3*32 + 9])) << 9;
            f[3] |= (*(image_center + i_y_a[3*32 + 10]*stride_image + i_x_a[3*32 + 10])
                        > *(image_center + i_y_b[3*32 + 10]*stride_image + i_x_b[3*32 + 10])) << 10;
            f[3] |= (*(image_center + i_y_a[3*32 + 11]*stride_image + i_x_a[3*32 + 11])
                        > *(image_center + i_y_b[3*32 + 11]*stride_image + i_x_b[3*32 + 11])) << 11;
            f[3] |= (*(image_center + i_y_a[3*32 + 12]*stride_image + i_x_a[3*32 + 12])
                        > *(image_center + i_y_b[3*32 + 12]*stride_image + i_x_b[3*32 + 12])) << 12;
            f[3] |= (*(image_center + i_y_a[3*32 + 13]*stride_image + i_x_a[3*32 + 13])
                        > *(image_center + i_y_b[3*32 + 13]*stride_image + i_x_b[3*32 + 13])) << 13;
            f[3] |= (*(image_center + i_y_a[3*32 + 14]*stride_image + i_x_a[3*32 + 14])
                        > *(image_center + i_y_b[3*32 + 14]*stride_image + i_x_b[3*32 + 14])) << 14;
            f[3] |= (*(image_center + i_y_a[3*32 + 15]*stride_image + i_x_a[3*32 + 15])
                        > *(image_center + i_y_b[3*32 + 15]*stride_image + i_x_b[3*32 + 15])) << 15;
            f[3] |= (*(image_center + i_y_a[3*32 + 16]*stride_image + i_x_a[3*32 + 16])
                        > *(image_center + i_y_b[3*32 + 16]*stride_image + i_x_b[3*32 + 16])) << 16;
            f[3] |= (*(image_center + i_y_a[3*32 + 17]*stride_image + i_x_a[3*32 + 17])
                        > *(image_center + i_y_b[3*32 + 17]*stride_image + i_x_b[3*32 + 17])) << 17;
            f[3] |= (*(image_center + i_y_a[3*32 + 18]*stride_image + i_x_a[3*32 + 18])
                        > *(image_center + i_y_b[3*32 + 18]*stride_image + i_x_b[3*32 + 18])) << 18;
            f[3] |= (*(image_center + i_y_a[3*32 + 19]*stride_image + i_x_a[3*32 + 19])
                        > *(image_center + i_y_b[3*32 + 19]*stride_image + i_x_b[3*32 + 19])) << 19;
            f[3] |= (*(image_center + i_y_a[3*32 + 20]*stride_image + i_x_a[3*32 + 20])
                        > *(image_center + i_y_b[3*32 + 20]*stride_image + i_x_b[3*32 + 20])) << 20;
            f[3] |= (*(image_center + i_y_a[3*32 + 21]*stride_image + i_x_a[3*32 + 21])
                        > *(image_center + i_y_b[3*32 + 21]*stride_image + i_x_b[3*32 + 21])) << 21;
            f[3] |= (*(image_center + i_y_a[3*32 + 22]*stride_image + i_x_a[3*32 + 22])
                        > *(image_center + i_y_b[3*32 + 22]*stride_image + i_x_b[3*32 + 22])) << 22;
            f[3] |= (*(image_center + i_y_a[3*32 + 23]*stride_image + i_x_a[3*32 + 23])
                        > *(image_center + i_y_b[3*32 + 23]*stride_image + i_x_b[3*32 + 23])) << 23;
            f[3] |= (*(image_center + i_y_a[3*32 + 24]*stride_image + i_x_a[3*32 + 24])
                        > *(image_center + i_y_b[3*32 + 24]*stride_image + i_x_b[3*32 + 24])) << 24;
            f[3] |= (*(image_center + i_y_a[3*32 + 25]*stride_image + i_x_a[3*32 + 25])
                        > *(image_center + i_y_b[3*32 + 25]*stride_image + i_x_b[3*32 + 25])) << 25;
            f[3] |= (*(image_center + i_y_a[3*32 + 26]*stride_image + i_x_a[3*32 + 26])
                        > *(image_center + i_y_b[3*32 + 26]*stride_image + i_x_b[3*32 + 26])) << 26;
            f[3] |= (*(image_center + i_y_a[3*32 + 27]*stride_image + i_x_a[3*32 + 27])
                        > *(image_center + i_y_b[3*32 + 27]*stride_image + i_x_b[3*32 + 27])) << 27;
            f[3] |= (*(image_center + i_y_a[3*32 + 28]*stride_image + i_x_a[3*32 + 28])
                        > *(image_center + i_y_b[3*32 + 28]*stride_image + i_x_b[3*32 + 28])) << 28;
            f[3] |= (*(image_center + i_y_a[3*32 + 29]*stride_image + i_x_a[3*32 + 29])
                        > *(image_center + i_y_b[3*32 + 29]*stride_image + i_x_b[3*32 + 29])) << 29;
            f[3] |= (*(image_center + i_y_a[3*32 + 30]*stride_image + i_x_a[3*32 + 30])
                        > *(image_center + i_y_b[3*32 + 30]*stride_image + i_x_b[3*32 + 30])) << 30;
            f[3] |= (*(image_center + i_y_a[3*32 + 31]*stride_image + i_x_a[3*32 + 31])
                        > *(image_center + i_y_b[3*32 + 31]*stride_image + i_x_b[3*32 + 31])) << 31;
            
            f[4] |= (*(image_center + i_y_a[4*32 + 0]*stride_image + i_x_a[4*32 + 0])
                        > *(image_center + i_y_b[4*32 + 0]*stride_image + i_x_b[4*32 + 0])) << 0;
            f[4] |= (*(image_center + i_y_a[4*32 + 1]*stride_image + i_x_a[4*32 + 1])
                        > *(image_center + i_y_b[4*32 + 1]*stride_image + i_x_b[4*32 + 1])) << 1;
            f[4] |= (*(image_center + i_y_a[4*32 + 2]*stride_image + i_x_a[4*32 + 2])
                        > *(image_center + i_y_b[4*32 + 2]*stride_image + i_x_b[4*32 + 2])) << 2;
            f[4] |= (*(image_center + i_y_a[4*32 + 3]*stride_image + i_x_a[4*32 + 3])
                        > *(image_center + i_y_b[4*32 + 3]*stride_image + i_x_b[4*32 + 3])) << 3;
            f[4] |= (*(image_center + i_y_a[4*32 + 4]*stride_image + i_x_a[4*32 + 4])
                        > *(image_center + i_y_b[4*32 + 4]*stride_image + i_x_b[4*32 + 4])) << 4;
            f[4] |= (*(image_center + i_y_a[4*32 + 5]*stride_image + i_x_a[4*32 + 5])
                        > *(image_center + i_y_b[4*32 + 5]*stride_image + i_x_b[4*32 + 5])) << 5;
            f[4] |= (*(image_center + i_y_a[4*32 + 6]*stride_image + i_x_a[4*32 + 6])
                        > *(image_center + i_y_b[4*32 + 6]*stride_image + i_x_b[4*32 + 6])) << 6;
            f[4] |= (*(image_center + i_y_a[4*32 + 7]*stride_image + i_x_a[4*32 + 7])
                        > *(image_center + i_y_b[4*32 + 7]*stride_image + i_x_b[4*32 + 7])) << 7;
            f[4] |= (*(image_center + i_y_a[4*32 + 8]*stride_image + i_x_a[4*32 + 8])
                        > *(image_center + i_y_b[4*32 + 8]*stride_image + i_x_b[4*32 + 8])) << 8;
            f[4] |= (*(image_center + i_y_a[4*32 + 9]*stride_image + i_x_a[4*32 + 9])
                        > *(image_center + i_y_b[4*32 + 9]*stride_image + i_x_b[4*32 + 9])) << 9;
            f[4] |= (*(image_center + i_y_a[4*32 + 10]*stride_image + i_x_a[4*32 + 10])
                        > *(image_center + i_y_b[4*32 + 10]*stride_image + i_x_b[4*32 + 10])) << 10;
            f[4] |= (*(image_center + i_y_a[4*32 + 11]*stride_image + i_x_a[4*32 + 11])
                        > *(image_center + i_y_b[4*32 + 11]*stride_image + i_x_b[4*32 + 11])) << 11;
            f[4] |= (*(image_center + i_y_a[4*32 + 12]*stride_image + i_x_a[4*32 + 12])
                        > *(image_center + i_y_b[4*32 + 12]*stride_image + i_x_b[4*32 + 12])) << 12;
            f[4] |= (*(image_center + i_y_a[4*32 + 13]*stride_image + i_x_a[4*32 + 13])
                        > *(image_center + i_y_b[4*32 + 13]*stride_image + i_x_b[4*32 + 13])) << 13;
            f[4] |= (*(image_center + i_y_a[4*32 + 14]*stride_image + i_x_a[4*32 + 14])
                        > *(image_center + i_y_b[4*32 + 14]*stride_image + i_x_b[4*32 + 14])) << 14;
            f[4] |= (*(image_center + i_y_a[4*32 + 15]*stride_image + i_x_a[4*32 + 15])
                        > *(image_center + i_y_b[4*32 + 15]*stride_image + i_x_b[4*32 + 15])) << 15;
            f[4] |= (*(image_center + i_y_a[4*32 + 16]*stride_image + i_x_a[4*32 + 16])
                        > *(image_center + i_y_b[4*32 + 16]*stride_image + i_x_b[4*32 + 16])) << 16;
            f[4] |= (*(image_center + i_y_a[4*32 + 17]*stride_image + i_x_a[4*32 + 17])
                        > *(image_center + i_y_b[4*32 + 17]*stride_image + i_x_b[4*32 + 17])) << 17;
            f[4] |= (*(image_center + i_y_a[4*32 + 18]*stride_image + i_x_a[4*32 + 18])
                        > *(image_center + i_y_b[4*32 + 18]*stride_image + i_x_b[4*32 + 18])) << 18;
            f[4] |= (*(image_center + i_y_a[4*32 + 19]*stride_image + i_x_a[4*32 + 19])
                        > *(image_center + i_y_b[4*32 + 19]*stride_image + i_x_b[4*32 + 19])) << 19;
            f[4] |= (*(image_center + i_y_a[4*32 + 20]*stride_image + i_x_a[4*32 + 20])
                        > *(image_center + i_y_b[4*32 + 20]*stride_image + i_x_b[4*32 + 20])) << 20;
            f[4] |= (*(image_center + i_y_a[4*32 + 21]*stride_image + i_x_a[4*32 + 21])
                        > *(image_center + i_y_b[4*32 + 21]*stride_image + i_x_b[4*32 + 21])) << 21;
            f[4] |= (*(image_center + i_y_a[4*32 + 22]*stride_image + i_x_a[4*32 + 22])
                        > *(image_center + i_y_b[4*32 + 22]*stride_image + i_x_b[4*32 + 22])) << 22;
            f[4] |= (*(image_center + i_y_a[4*32 + 23]*stride_image + i_x_a[4*32 + 23])
                        > *(image_center + i_y_b[4*32 + 23]*stride_image + i_x_b[4*32 + 23])) << 23;
            f[4] |= (*(image_center + i_y_a[4*32 + 24]*stride_image + i_x_a[4*32 + 24])
                        > *(image_center + i_y_b[4*32 + 24]*stride_image + i_x_b[4*32 + 24])) << 24;
            f[4] |= (*(image_center + i_y_a[4*32 + 25]*stride_image + i_x_a[4*32 + 25])
                        > *(image_center + i_y_b[4*32 + 25]*stride_image + i_x_b[4*32 + 25])) << 25;
            f[4] |= (*(image_center + i_y_a[4*32 + 26]*stride_image + i_x_a[4*32 + 26])
                        > *(image_center + i_y_b[4*32 + 26]*stride_image + i_x_b[4*32 + 26])) << 26;
            f[4] |= (*(image_center + i_y_a[4*32 + 27]*stride_image + i_x_a[4*32 + 27])
                        > *(image_center + i_y_b[4*32 + 27]*stride_image + i_x_b[4*32 + 27])) << 27;
            f[4] |= (*(image_center + i_y_a[4*32 + 28]*stride_image + i_x_a[4*32 + 28])
                        > *(image_center + i_y_b[4*32 + 28]*stride_image + i_x_b[4*32 + 28])) << 28;
            f[4] |= (*(image_center + i_y_a[4*32 + 29]*stride_image + i_x_a[4*32 + 29])
                        > *(image_center + i_y_b[4*32 + 29]*stride_image + i_x_b[4*32 + 29])) << 29;
            f[4] |= (*(image_center + i_y_a[4*32 + 30]*stride_image + i_x_a[4*32 + 30])
                        > *(image_center + i_y_b[4*32 + 30]*stride_image + i_x_b[4*32 + 30])) << 30;
            f[4] |= (*(image_center + i_y_a[4*32 + 31]*stride_image + i_x_a[4*32 + 31])
                        > *(image_center + i_y_b[4*32 + 31]*stride_image + i_x_b[4*32 + 31])) << 31;
            
            f[5] |= (*(image_center + i_y_a[5*32 + 0]*stride_image + i_x_a[5*32 + 0])
                        > *(image_center + i_y_b[5*32 + 0]*stride_image + i_x_b[5*32 + 0])) << 0;
            f[5] |= (*(image_center + i_y_a[5*32 + 1]*stride_image + i_x_a[5*32 + 1])
                        > *(image_center + i_y_b[5*32 + 1]*stride_image + i_x_b[5*32 + 1])) << 1;
            f[5] |= (*(image_center + i_y_a[5*32 + 2]*stride_image + i_x_a[5*32 + 2])
                        > *(image_center + i_y_b[5*32 + 2]*stride_image + i_x_b[5*32 + 2])) << 2;
            f[5] |= (*(image_center + i_y_a[5*32 + 3]*stride_image + i_x_a[5*32 + 3])
                        > *(image_center + i_y_b[5*32 + 3]*stride_image + i_x_b[5*32 + 3])) << 3;
            f[5] |= (*(image_center + i_y_a[5*32 + 4]*stride_image + i_x_a[5*32 + 4])
                        > *(image_center + i_y_b[5*32 + 4]*stride_image + i_x_b[5*32 + 4])) << 4;
            f[5] |= (*(image_center + i_y_a[5*32 + 5]*stride_image + i_x_a[5*32 + 5])
                        > *(image_center + i_y_b[5*32 + 5]*stride_image + i_x_b[5*32 + 5])) << 5;
            f[5] |= (*(image_center + i_y_a[5*32 + 6]*stride_image + i_x_a[5*32 + 6])
                        > *(image_center + i_y_b[5*32 + 6]*stride_image + i_x_b[5*32 + 6])) << 6;
            f[5] |= (*(image_center + i_y_a[5*32 + 7]*stride_image + i_x_a[5*32 + 7])
                        > *(image_center + i_y_b[5*32 + 7]*stride_image + i_x_b[5*32 + 7])) << 7;
            f[5] |= (*(image_center + i_y_a[5*32 + 8]*stride_image + i_x_a[5*32 + 8])
                        > *(image_center + i_y_b[5*32 + 8]*stride_image + i_x_b[5*32 + 8])) << 8;
            f[5] |= (*(image_center + i_y_a[5*32 + 9]*stride_image + i_x_a[5*32 + 9])
                        > *(image_center + i_y_b[5*32 + 9]*stride_image + i_x_b[5*32 + 9])) << 9;
            f[5] |= (*(image_center + i_y_a[5*32 + 10]*stride_image + i_x_a[5*32 + 10])
                        > *(image_center + i_y_b[5*32 + 10]*stride_image + i_x_b[5*32 + 10])) << 10;
            f[5] |= (*(image_center + i_y_a[5*32 + 11]*stride_image + i_x_a[5*32 + 11])
                        > *(image_center + i_y_b[5*32 + 11]*stride_image + i_x_b[5*32 + 11])) << 11;
            f[5] |= (*(image_center + i_y_a[5*32 + 12]*stride_image + i_x_a[5*32 + 12])
                        > *(image_center + i_y_b[5*32 + 12]*stride_image + i_x_b[5*32 + 12])) << 12;
            f[5] |= (*(image_center + i_y_a[5*32 + 13]*stride_image + i_x_a[5*32 + 13])
                        > *(image_center + i_y_b[5*32 + 13]*stride_image + i_x_b[5*32 + 13])) << 13;
            f[5] |= (*(image_center + i_y_a[5*32 + 14]*stride_image + i_x_a[5*32 + 14])
                        > *(image_center + i_y_b[5*32 + 14]*stride_image + i_x_b[5*32 + 14])) << 14;
            f[5] |= (*(image_center + i_y_a[5*32 + 15]*stride_image + i_x_a[5*32 + 15])
                        > *(image_center + i_y_b[5*32 + 15]*stride_image + i_x_b[5*32 + 15])) << 15;
            f[5] |= (*(image_center + i_y_a[5*32 + 16]*stride_image + i_x_a[5*32 + 16])
                        > *(image_center + i_y_b[5*32 + 16]*stride_image + i_x_b[5*32 + 16])) << 16;
            f[5] |= (*(image_center + i_y_a[5*32 + 17]*stride_image + i_x_a[5*32 + 17])
                        > *(image_center + i_y_b[5*32 + 17]*stride_image + i_x_b[5*32 + 17])) << 17;
            f[5] |= (*(image_center + i_y_a[5*32 + 18]*stride_image + i_x_a[5*32 + 18])
                        > *(image_center + i_y_b[5*32 + 18]*stride_image + i_x_b[5*32 + 18])) << 18;
            f[5] |= (*(image_center + i_y_a[5*32 + 19]*stride_image + i_x_a[5*32 + 19])
                        > *(image_center + i_y_b[5*32 + 19]*stride_image + i_x_b[5*32 + 19])) << 19;
            f[5] |= (*(image_center + i_y_a[5*32 + 20]*stride_image + i_x_a[5*32 + 20])
                        > *(image_center + i_y_b[5*32 + 20]*stride_image + i_x_b[5*32 + 20])) << 20;
            f[5] |= (*(image_center + i_y_a[5*32 + 21]*stride_image + i_x_a[5*32 + 21])
                        > *(image_center + i_y_b[5*32 + 21]*stride_image + i_x_b[5*32 + 21])) << 21;
            f[5] |= (*(image_center + i_y_a[5*32 + 22]*stride_image + i_x_a[5*32 + 22])
                        > *(image_center + i_y_b[5*32 + 22]*stride_image + i_x_b[5*32 + 22])) << 22;
            f[5] |= (*(image_center + i_y_a[5*32 + 23]*stride_image + i_x_a[5*32 + 23])
                        > *(image_center + i_y_b[5*32 + 23]*stride_image + i_x_b[5*32 + 23])) << 23;
            f[5] |= (*(image_center + i_y_a[5*32 + 24]*stride_image + i_x_a[5*32 + 24])
                        > *(image_center + i_y_b[5*32 + 24]*stride_image + i_x_b[5*32 + 24])) << 24;
            f[5] |= (*(image_center + i_y_a[5*32 + 25]*stride_image + i_x_a[5*32 + 25])
                        > *(image_center + i_y_b[5*32 + 25]*stride_image + i_x_b[5*32 + 25])) << 25;
            f[5] |= (*(image_center + i_y_a[5*32 + 26]*stride_image + i_x_a[5*32 + 26])
                        > *(image_center + i_y_b[5*32 + 26]*stride_image + i_x_b[5*32 + 26])) << 26;
            f[5] |= (*(image_center + i_y_a[5*32 + 27]*stride_image + i_x_a[5*32 + 27])
                        > *(image_center + i_y_b[5*32 + 27]*stride_image + i_x_b[5*32 + 27])) << 27;
            f[5] |= (*(image_center + i_y_a[5*32 + 28]*stride_image + i_x_a[5*32 + 28])
                        > *(image_center + i_y_b[5*32 + 28]*stride_image + i_x_b[5*32 + 28])) << 28;
            f[5] |= (*(image_center + i_y_a[5*32 + 29]*stride_image + i_x_a[5*32 + 29])
                        > *(image_center + i_y_b[5*32 + 29]*stride_image + i_x_b[5*32 + 29])) << 29;
            f[5] |= (*(image_center + i_y_a[5*32 + 30]*stride_image + i_x_a[5*32 + 30])
                        > *(image_center + i_y_b[5*32 + 30]*stride_image + i_x_b[5*32 + 30])) << 30;
            f[5] |= (*(image_center + i_y_a[5*32 + 31]*stride_image + i_x_a[5*32 + 31])
                        > *(image_center + i_y_b[5*32 + 31]*stride_image + i_x_b[5*32 + 31])) << 31;
            
            f[6] |= (*(image_center + i_y_a[6*32 + 0]*stride_image + i_x_a[6*32 + 0])
                        > *(image_center + i_y_b[6*32 + 0]*stride_image + i_x_b[6*32 + 0])) << 0;
            f[6] |= (*(image_center + i_y_a[6*32 + 1]*stride_image + i_x_a[6*32 + 1])
                        > *(image_center + i_y_b[6*32 + 1]*stride_image + i_x_b[6*32 + 1])) << 1;
            f[6] |= (*(image_center + i_y_a[6*32 + 2]*stride_image + i_x_a[6*32 + 2])
                        > *(image_center + i_y_b[6*32 + 2]*stride_image + i_x_b[6*32 + 2])) << 2;
            f[6] |= (*(image_center + i_y_a[6*32 + 3]*stride_image + i_x_a[6*32 + 3])
                        > *(image_center + i_y_b[6*32 + 3]*stride_image + i_x_b[6*32 + 3])) << 3;
            f[6] |= (*(image_center + i_y_a[6*32 + 4]*stride_image + i_x_a[6*32 + 4])
                        > *(image_center + i_y_b[6*32 + 4]*stride_image + i_x_b[6*32 + 4])) << 4;
            f[6] |= (*(image_center + i_y_a[6*32 + 5]*stride_image + i_x_a[6*32 + 5])
                        > *(image_center + i_y_b[6*32 + 5]*stride_image + i_x_b[6*32 + 5])) << 5;
            f[6] |= (*(image_center + i_y_a[6*32 + 6]*stride_image + i_x_a[6*32 + 6])
                        > *(image_center + i_y_b[6*32 + 6]*stride_image + i_x_b[6*32 + 6])) << 6;
            f[6] |= (*(image_center + i_y_a[6*32 + 7]*stride_image + i_x_a[6*32 + 7])
                        > *(image_center + i_y_b[6*32 + 7]*stride_image + i_x_b[6*32 + 7])) << 7;
            f[6] |= (*(image_center + i_y_a[6*32 + 8]*stride_image + i_x_a[6*32 + 8])
                        > *(image_center + i_y_b[6*32 + 8]*stride_image + i_x_b[6*32 + 8])) << 8;
            f[6] |= (*(image_center + i_y_a[6*32 + 9]*stride_image + i_x_a[6*32 + 9])
                        > *(image_center + i_y_b[6*32 + 9]*stride_image + i_x_b[6*32 + 9])) << 9;
            f[6] |= (*(image_center + i_y_a[6*32 + 10]*stride_image + i_x_a[6*32 + 10])
                        > *(image_center + i_y_b[6*32 + 10]*stride_image + i_x_b[6*32 + 10])) << 10;
            f[6] |= (*(image_center + i_y_a[6*32 + 11]*stride_image + i_x_a[6*32 + 11])
                        > *(image_center + i_y_b[6*32 + 11]*stride_image + i_x_b[6*32 + 11])) << 11;
            f[6] |= (*(image_center + i_y_a[6*32 + 12]*stride_image + i_x_a[6*32 + 12])
                        > *(image_center + i_y_b[6*32 + 12]*stride_image + i_x_b[6*32 + 12])) << 12;
            f[6] |= (*(image_center + i_y_a[6*32 + 13]*stride_image + i_x_a[6*32 + 13])
                        > *(image_center + i_y_b[6*32 + 13]*stride_image + i_x_b[6*32 + 13])) << 13;
            f[6] |= (*(image_center + i_y_a[6*32 + 14]*stride_image + i_x_a[6*32 + 14])
                        > *(image_center + i_y_b[6*32 + 14]*stride_image + i_x_b[6*32 + 14])) << 14;
            f[6] |= (*(image_center + i_y_a[6*32 + 15]*stride_image + i_x_a[6*32 + 15])
                        > *(image_center + i_y_b[6*32 + 15]*stride_image + i_x_b[6*32 + 15])) << 15;
            f[6] |= (*(image_center + i_y_a[6*32 + 16]*stride_image + i_x_a[6*32 + 16])
                        > *(image_center + i_y_b[6*32 + 16]*stride_image + i_x_b[6*32 + 16])) << 16;
            f[6] |= (*(image_center + i_y_a[6*32 + 17]*stride_image + i_x_a[6*32 + 17])
                        > *(image_center + i_y_b[6*32 + 17]*stride_image + i_x_b[6*32 + 17])) << 17;
            f[6] |= (*(image_center + i_y_a[6*32 + 18]*stride_image + i_x_a[6*32 + 18])
                        > *(image_center + i_y_b[6*32 + 18]*stride_image + i_x_b[6*32 + 18])) << 18;
            f[6] |= (*(image_center + i_y_a[6*32 + 19]*stride_image + i_x_a[6*32 + 19])
                        > *(image_center + i_y_b[6*32 + 19]*stride_image + i_x_b[6*32 + 19])) << 19;
            f[6] |= (*(image_center + i_y_a[6*32 + 20]*stride_image + i_x_a[6*32 + 20])
                        > *(image_center + i_y_b[6*32 + 20]*stride_image + i_x_b[6*32 + 20])) << 20;
            f[6] |= (*(image_center + i_y_a[6*32 + 21]*stride_image + i_x_a[6*32 + 21])
                        > *(image_center + i_y_b[6*32 + 21]*stride_image + i_x_b[6*32 + 21])) << 21;
            f[6] |= (*(image_center + i_y_a[6*32 + 22]*stride_image + i_x_a[6*32 + 22])
                        > *(image_center + i_y_b[6*32 + 22]*stride_image + i_x_b[6*32 + 22])) << 22;
            f[6] |= (*(image_center + i_y_a[6*32 + 23]*stride_image + i_x_a[6*32 + 23])
                        > *(image_center + i_y_b[6*32 + 23]*stride_image + i_x_b[6*32 + 23])) << 23;
            f[6] |= (*(image_center + i_y_a[6*32 + 24]*stride_image + i_x_a[6*32 + 24])
                        > *(image_center + i_y_b[6*32 + 24]*stride_image + i_x_b[6*32 + 24])) << 24;
            f[6] |= (*(image_center + i_y_a[6*32 + 25]*stride_image + i_x_a[6*32 + 25])
                        > *(image_center + i_y_b[6*32 + 25]*stride_image + i_x_b[6*32 + 25])) << 25;
            f[6] |= (*(image_center + i_y_a[6*32 + 26]*stride_image + i_x_a[6*32 + 26])
                        > *(image_center + i_y_b[6*32 + 26]*stride_image + i_x_b[6*32 + 26])) << 26;
            f[6] |= (*(image_center + i_y_a[6*32 + 27]*stride_image + i_x_a[6*32 + 27])
                        > *(image_center + i_y_b[6*32 + 27]*stride_image + i_x_b[6*32 + 27])) << 27;
            f[6] |= (*(image_center + i_y_a[6*32 + 28]*stride_image + i_x_a[6*32 + 28])
                        > *(image_center + i_y_b[6*32 + 28]*stride_image + i_x_b[6*32 + 28])) << 28;
            f[6] |= (*(image_center + i_y_a[6*32 + 29]*stride_image + i_x_a[6*32 + 29])
                        > *(image_center + i_y_b[6*32 + 29]*stride_image + i_x_b[6*32 + 29])) << 29;
            f[6] |= (*(image_center + i_y_a[6*32 + 30]*stride_image + i_x_a[6*32 + 30])
                        > *(image_center + i_y_b[6*32 + 30]*stride_image + i_x_b[6*32 + 30])) << 30;
            f[6] |= (*(image_center + i_y_a[6*32 + 31]*stride_image + i_x_a[6*32 + 31])
                        > *(image_center + i_y_b[6*32 + 31]*stride_image + i_x_b[6*32 + 31])) << 31;
            
            f[7] |= (*(image_center + i_y_a[7*32 + 0]*stride_image + i_x_a[7*32 + 0])
                        > *(image_center + i_y_b[7*32 + 0]*stride_image + i_x_b[7*32 + 0])) << 0;
            f[7] |= (*(image_center + i_y_a[7*32 + 1]*stride_image + i_x_a[7*32 + 1])
                        > *(image_center + i_y_b[7*32 + 1]*stride_image + i_x_b[7*32 + 1])) << 1;
            f[7] |= (*(image_center + i_y_a[7*32 + 2]*stride_image + i_x_a[7*32 + 2])
                        > *(image_center + i_y_b[7*32 + 2]*stride_image + i_x_b[7*32 + 2])) << 2;
            f[7] |= (*(image_center + i_y_a[7*32 + 3]*stride_image + i_x_a[7*32 + 3])
                        > *(image_center + i_y_b[7*32 + 3]*stride_image + i_x_b[7*32 + 3])) << 3;
            f[7] |= (*(image_center + i_y_a[7*32 + 4]*stride_image + i_x_a[7*32 + 4])
                        > *(image_center + i_y_b[7*32 + 4]*stride_image + i_x_b[7*32 + 4])) << 4;
            f[7] |= (*(image_center + i_y_a[7*32 + 5]*stride_image + i_x_a[7*32 + 5])
                        > *(image_center + i_y_b[7*32 + 5]*stride_image + i_x_b[7*32 + 5])) << 5;
            f[7] |= (*(image_center + i_y_a[7*32 + 6]*stride_image + i_x_a[7*32 + 6])
                        > *(image_center + i_y_b[7*32 + 6]*stride_image + i_x_b[7*32 + 6])) << 6;
            f[7] |= (*(image_center + i_y_a[7*32 + 7]*stride_image + i_x_a[7*32 + 7])
                        > *(image_center + i_y_b[7*32 + 7]*stride_image + i_x_b[7*32 + 7])) << 7;
            f[7] |= (*(image_center + i_y_a[7*32 + 8]*stride_image + i_x_a[7*32 + 8])
                        > *(image_center + i_y_b[7*32 + 8]*stride_image + i_x_b[7*32 + 8])) << 8;
            f[7] |= (*(image_center + i_y_a[7*32 + 9]*stride_image + i_x_a[7*32 + 9])
                        > *(image_center + i_y_b[7*32 + 9]*stride_image + i_x_b[7*32 + 9])) << 9;
            f[7] |= (*(image_center + i_y_a[7*32 + 10]*stride_image + i_x_a[7*32 + 10])
                        > *(image_center + i_y_b[7*32 + 10]*stride_image + i_x_b[7*32 + 10])) << 10;
            f[7] |= (*(image_center + i_y_a[7*32 + 11]*stride_image + i_x_a[7*32 + 11])
                        > *(image_center + i_y_b[7*32 + 11]*stride_image + i_x_b[7*32 + 11])) << 11;
            f[7] |= (*(image_center + i_y_a[7*32 + 12]*stride_image + i_x_a[7*32 + 12])
                        > *(image_center + i_y_b[7*32 + 12]*stride_image + i_x_b[7*32 + 12])) << 12;
            f[7] |= (*(image_center + i_y_a[7*32 + 13]*stride_image + i_x_a[7*32 + 13])
                        > *(image_center + i_y_b[7*32 + 13]*stride_image + i_x_b[7*32 + 13])) << 13;
            f[7] |= (*(image_center + i_y_a[7*32 + 14]*stride_image + i_x_a[7*32 + 14])
                        > *(image_center + i_y_b[7*32 + 14]*stride_image + i_x_b[7*32 + 14])) << 14;
            f[7] |= (*(image_center + i_y_a[7*32 + 15]*stride_image + i_x_a[7*32 + 15])
                        > *(image_center + i_y_b[7*32 + 15]*stride_image + i_x_b[7*32 + 15])) << 15;
            f[7] |= (*(image_center + i_y_a[7*32 + 16]*stride_image + i_x_a[7*32 + 16])
                        > *(image_center + i_y_b[7*32 + 16]*stride_image + i_x_b[7*32 + 16])) << 16;
            f[7] |= (*(image_center + i_y_a[7*32 + 17]*stride_image + i_x_a[7*32 + 17])
                        > *(image_center + i_y_b[7*32 + 17]*stride_image + i_x_b[7*32 + 17])) << 17;
            f[7] |= (*(image_center + i_y_a[7*32 + 18]*stride_image + i_x_a[7*32 + 18])
                        > *(image_center + i_y_b[7*32 + 18]*stride_image + i_x_b[7*32 + 18])) << 18;
            f[7] |= (*(image_center + i_y_a[7*32 + 19]*stride_image + i_x_a[7*32 + 19])
                        > *(image_center + i_y_b[7*32 + 19]*stride_image + i_x_b[7*32 + 19])) << 19;
            f[7] |= (*(image_center + i_y_a[7*32 + 20]*stride_image + i_x_a[7*32 + 20])
                        > *(image_center + i_y_b[7*32 + 20]*stride_image + i_x_b[7*32 + 20])) << 20;
            f[7] |= (*(image_center + i_y_a[7*32 + 21]*stride_image + i_x_a[7*32 + 21])
                        > *(image_center + i_y_b[7*32 + 21]*stride_image + i_x_b[7*32 + 21])) << 21;
            f[7] |= (*(image_center + i_y_a[7*32 + 22]*stride_image + i_x_a[7*32 + 22])
                        > *(image_center + i_y_b[7*32 + 22]*stride_image + i_x_b[7*32 + 22])) << 22;
            f[7] |= (*(image_center + i_y_a[7*32 + 23]*stride_image + i_x_a[7*32 + 23])
                        > *(image_center + i_y_b[7*32 + 23]*stride_image + i_x_b[7*32 + 23])) << 23;
            f[7] |= (*(image_center + i_y_a[7*32 + 24]*stride_image + i_x_a[7*32 + 24])
                        > *(image_center + i_y_b[7*32 + 24]*stride_image + i_x_b[7*32 + 24])) << 24;
            f[7] |= (*(image_center + i_y_a[7*32 + 25]*stride_image + i_x_a[7*32 + 25])
                        > *(image_center + i_y_b[7*32 + 25]*stride_image + i_x_b[7*32 + 25])) << 25;
            f[7] |= (*(image_center + i_y_a[7*32 + 26]*stride_image + i_x_a[7*32 + 26])
                        > *(image_center + i_y_b[7*32 + 26]*stride_image + i_x_b[7*32 + 26])) << 26;
            f[7] |= (*(image_center + i_y_a[7*32 + 27]*stride_image + i_x_a[7*32 + 27])
                        > *(image_center + i_y_b[7*32 + 27]*stride_image + i_x_b[7*32 + 27])) << 27;
            f[7] |= (*(image_center + i_y_a[7*32 + 28]*stride_image + i_x_a[7*32 + 28])
                        > *(image_center + i_y_b[7*32 + 28]*stride_image + i_x_b[7*32 + 28])) << 28;
            f[7] |= (*(image_center + i_y_a[7*32 + 29]*stride_image + i_x_a[7*32 + 29])
                        > *(image_center + i_y_b[7*32 + 29]*stride_image + i_x_b[7*32 + 29])) << 29;
            f[7] |= (*(image_center + i_y_a[7*32 + 30]*stride_image + i_x_a[7*32 + 30])
                        > *(image_center + i_y_b[7*32 + 30]*stride_image + i_x_b[7*32 + 30])) << 30;
            f[7] |= (*(image_center + i_y_a[7*32 + 31]*stride_image + i_x_a[7*32 + 31])
                        > *(image_center + i_y_b[7*32 + 31]*stride_image + i_x_b[7*32 + 31])) << 31;
            
            f[8] |= (*(image_center + i_y_a[8*32 + 0]*stride_image + i_x_a[8*32 + 0])
                        > *(image_center + i_y_b[8*32 + 0]*stride_image + i_x_b[8*32 + 0])) << 0;
            f[8] |= (*(image_center + i_y_a[8*32 + 1]*stride_image + i_x_a[8*32 + 1])
                        > *(image_center + i_y_b[8*32 + 1]*stride_image + i_x_b[8*32 + 1])) << 1;
            f[8] |= (*(image_center + i_y_a[8*32 + 2]*stride_image + i_x_a[8*32 + 2])
                        > *(image_center + i_y_b[8*32 + 2]*stride_image + i_x_b[8*32 + 2])) << 2;
            f[8] |= (*(image_center + i_y_a[8*32 + 3]*stride_image + i_x_a[8*32 + 3])
                        > *(image_center + i_y_b[8*32 + 3]*stride_image + i_x_b[8*32 + 3])) << 3;
            f[8] |= (*(image_center + i_y_a[8*32 + 4]*stride_image + i_x_a[8*32 + 4])
                        > *(image_center + i_y_b[8*32 + 4]*stride_image + i_x_b[8*32 + 4])) << 4;
            f[8] |= (*(image_center + i_y_a[8*32 + 5]*stride_image + i_x_a[8*32 + 5])
                        > *(image_center + i_y_b[8*32 + 5]*stride_image + i_x_b[8*32 + 5])) << 5;
            f[8] |= (*(image_center + i_y_a[8*32 + 6]*stride_image + i_x_a[8*32 + 6])
                        > *(image_center + i_y_b[8*32 + 6]*stride_image + i_x_b[8*32 + 6])) << 6;
            f[8] |= (*(image_center + i_y_a[8*32 + 7]*stride_image + i_x_a[8*32 + 7])
                        > *(image_center + i_y_b[8*32 + 7]*stride_image + i_x_b[8*32 + 7])) << 7;
            f[8] |= (*(image_center + i_y_a[8*32 + 8]*stride_image + i_x_a[8*32 + 8])
                        > *(image_center + i_y_b[8*32 + 8]*stride_image + i_x_b[8*32 + 8])) << 8;
            f[8] |= (*(image_center + i_y_a[8*32 + 9]*stride_image + i_x_a[8*32 + 9])
                        > *(image_center + i_y_b[8*32 + 9]*stride_image + i_x_b[8*32 + 9])) << 9;
            f[8] |= (*(image_center + i_y_a[8*32 + 10]*stride_image + i_x_a[8*32 + 10])
                        > *(image_center + i_y_b[8*32 + 10]*stride_image + i_x_b[8*32 + 10])) << 10;
            f[8] |= (*(image_center + i_y_a[8*32 + 11]*stride_image + i_x_a[8*32 + 11])
                        > *(image_center + i_y_b[8*32 + 11]*stride_image + i_x_b[8*32 + 11])) << 11;
            f[8] |= (*(image_center + i_y_a[8*32 + 12]*stride_image + i_x_a[8*32 + 12])
                        > *(image_center + i_y_b[8*32 + 12]*stride_image + i_x_b[8*32 + 12])) << 12;
            f[8] |= (*(image_center + i_y_a[8*32 + 13]*stride_image + i_x_a[8*32 + 13])
                        > *(image_center + i_y_b[8*32 + 13]*stride_image + i_x_b[8*32 + 13])) << 13;
            f[8] |= (*(image_center + i_y_a[8*32 + 14]*stride_image + i_x_a[8*32 + 14])
                        > *(image_center + i_y_b[8*32 + 14]*stride_image + i_x_b[8*32 + 14])) << 14;
            f[8] |= (*(image_center + i_y_a[8*32 + 15]*stride_image + i_x_a[8*32 + 15])
                        > *(image_center + i_y_b[8*32 + 15]*stride_image + i_x_b[8*32 + 15])) << 15;
            f[8] |= (*(image_center + i_y_a[8*32 + 16]*stride_image + i_x_a[8*32 + 16])
                        > *(image_center + i_y_b[8*32 + 16]*stride_image + i_x_b[8*32 + 16])) << 16;
            f[8] |= (*(image_center + i_y_a[8*32 + 17]*stride_image + i_x_a[8*32 + 17])
                        > *(image_center + i_y_b[8*32 + 17]*stride_image + i_x_b[8*32 + 17])) << 17;
            f[8] |= (*(image_center + i_y_a[8*32 + 18]*stride_image + i_x_a[8*32 + 18])
                        > *(image_center + i_y_b[8*32 + 18]*stride_image + i_x_b[8*32 + 18])) << 18;
            f[8] |= (*(image_center + i_y_a[8*32 + 19]*stride_image + i_x_a[8*32 + 19])
                        > *(image_center + i_y_b[8*32 + 19]*stride_image + i_x_b[8*32 + 19])) << 19;
            f[8] |= (*(image_center + i_y_a[8*32 + 20]*stride_image + i_x_a[8*32 + 20])
                        > *(image_center + i_y_b[8*32 + 20]*stride_image + i_x_b[8*32 + 20])) << 20;
            f[8] |= (*(image_center + i_y_a[8*32 + 21]*stride_image + i_x_a[8*32 + 21])
                        > *(image_center + i_y_b[8*32 + 21]*stride_image + i_x_b[8*32 + 21])) << 21;
            f[8] |= (*(image_center + i_y_a[8*32 + 22]*stride_image + i_x_a[8*32 + 22])
                        > *(image_center + i_y_b[8*32 + 22]*stride_image + i_x_b[8*32 + 22])) << 22;
            f[8] |= (*(image_center + i_y_a[8*32 + 23]*stride_image + i_x_a[8*32 + 23])
                        > *(image_center + i_y_b[8*32 + 23]*stride_image + i_x_b[8*32 + 23])) << 23;
            f[8] |= (*(image_center + i_y_a[8*32 + 24]*stride_image + i_x_a[8*32 + 24])
                        > *(image_center + i_y_b[8*32 + 24]*stride_image + i_x_b[8*32 + 24])) << 24;
            f[8] |= (*(image_center + i_y_a[8*32 + 25]*stride_image + i_x_a[8*32 + 25])
                        > *(image_center + i_y_b[8*32 + 25]*stride_image + i_x_b[8*32 + 25])) << 25;
            f[8] |= (*(image_center + i_y_a[8*32 + 26]*stride_image + i_x_a[8*32 + 26])
                        > *(image_center + i_y_b[8*32 + 26]*stride_image + i_x_b[8*32 + 26])) << 26;
            f[8] |= (*(image_center + i_y_a[8*32 + 27]*stride_image + i_x_a[8*32 + 27])
                        > *(image_center + i_y_b[8*32 + 27]*stride_image + i_x_b[8*32 + 27])) << 27;
            f[8] |= (*(image_center + i_y_a[8*32 + 28]*stride_image + i_x_a[8*32 + 28])
                        > *(image_center + i_y_b[8*32 + 28]*stride_image + i_x_b[8*32 + 28])) << 28;
            f[8] |= (*(image_center + i_y_a[8*32 + 29]*stride_image + i_x_a[8*32 + 29])
                        > *(image_center + i_y_b[8*32 + 29]*stride_image + i_x_b[8*32 + 29])) << 29;
            f[8] |= (*(image_center + i_y_a[8*32 + 30]*stride_image + i_x_a[8*32 + 30])
                        > *(image_center + i_y_b[8*32 + 30]*stride_image + i_x_b[8*32 + 30])) << 30;
            f[8] |= (*(image_center + i_y_a[8*32 + 31]*stride_image + i_x_a[8*32 + 31])
                        > *(image_center + i_y_b[8*32 + 31]*stride_image + i_x_b[8*32 + 31])) << 31;
            
            
            _mm_store_si128((__m128i*)(bd.memptr() + 0*n_rows),_mm_load_si128((const __m128i *)(f + 0*64)));
            _mm_store_si128((__m128i*)(bd.memptr() + 1*n_rows),_mm_load_si128((const __m128i *)(f + 1*64)));
            _mm_store_si128((__m128i*)(bd.memptr() + 2*n_rows),_mm_load_si128((const __m128i *)(f + 2*64)));
            _mm_store_si128((__m128i*)(bd.memptr() + 3*n_rows),_mm_load_si128((const __m128i *)(f + 3*64)));
            

        }
    }
}


int BRIEF::diag_length_pattern = 17;
int BRIEF::gaussian_bit_pattern_31_x_a[256] = { 8,4,-11,7,2,1,-2,-13,-13,10,-13,-11,7,-4,-13,-9,12,-3,-6,11,4,5,3,-8,-2,-13,-7,-4,-10,5,5,1,9,4,2,-4,-8,4,0,-13,-3,-6,8,0,7,-13,10,-6,10,-13,-13,3,5,-1,3,2,-13,-13,-13,-7,6,-9,-2,-12,3,-7,-3,2,-11,-1,5,-4,-9,-12,10,7,-7,-4,7,-7,-13,-3,7,-13,1,2,-4,-1,7,1,9,-1,-13,7,12,6,5,2,3,2,9,-8,-11,1,6,2,6,3,7,-11,-10,-5,-10,8,4,-10,4,-2,-5,7,-9,-5,8,-9,1,7,-2,11,-12,3,5,0,-9,0,-1,5,3,-13,-5,-4,6,-7,-13,1,4,-2,2,-2,4,-6,-3,7,4,-13,7,7,-7,-8,-13,2,10,-6,8,2,-11,-12,-11,5,-2,-1,-13,-10,-3,2,-9,-4,-4,-6,6,-13,11,7,-1,-4,-7,-13,-7,-8,-5,-13,1,1,9,5,-1,-9,-1,-13,8,2,7,-10,-10,4,3,-4,5,4,-9,0,-12,3,-10,8,-8,2,10,6,-7,-3,-1,-3,-8,4,2,6,3,11,-3,4,2,-10,-13,-13,6,0,-13,-9,-13,5,2,-1,9,11,3,-1,3,-13,5,8,7,-10,7,9,7,-1};
int BRIEF::gaussian_bit_pattern_31_y_a[256] = {-3,2,9,-12,-13,-7,-10,-13,-3,4,-8,7,7,-5,2,0,-6,6,-13,-13,7,-3,-7,-7,11,12,3,2,-12,-12,-6,0,11,7,-1,-12,-5,11,-8,-2,-2,9,12,9,-5,-6,7,-3,-9,8,0,3,7,7,-10,-4,0,-7,3,12,-10,-1,-5,5,-10,-7,-2,9,-13,6,-3,-13,-6,-10,2,12,-13,9,-1,6,11,7,-8,-7,-3,-6,3,-13,1,-1,1,-9,-13,7,-5,3,-13,-12,8,6,-12,4,12,12,-9,3,3,-3,8,-5,11,-8,5,-1,-6,12,-2,0,-8,-6,-13,-13,-8,-11,-8,-4,1,-6,-9,7,5,-4,12,7,2,11,5,-4,9,-7,5,6,6,-10,1,-2,-12,-13,1,-10,-13,5,-2,9,1,-8,-4,11,6,4,-5,-5,-3,-12,-2,-13,0,-3,-13,-8,-11,-2,9,-3,-13,6,12,-11,-3,11,11,-5,12,-8,1,-12,-2,5,-1,7,5,0,12,-8,11,-3,-10,1,-11,-13,-13,-10,-8,-6,12,2,-13,-13,9,3,1,2,-10,-13,-12,2,6,8,10,-9,-13,-7,-2,2,-5,-9,-1,-1,0,-11,-4,-6,7,12,0,-1,3,8,-6,-9,7,-6,5,-3,0,4,-6,0,8,9,-4,4,3,-7,0,-6};
int BRIEF::gaussian_bit_pattern_31_x_b[256] = {9,7,-8,12,2,1,-2,-11,-12,11,-8,-9,12,-3,-12,-7,12,-2,-4,12,5,10,6,-6,-1,-8,-5,-3,-6,6,7,4,11,4,4,-2,-7,9,1,-8,-2,-4,10,1,11,-11,12,-6,12,-8,-8,7,10,1,5,3,-13,-12,-11,-4,12,-7,0,-7,8,-4,-1,5,-5,0,5,-4,-9,-8,12,12,-6,-3,12,-5,-12,-2,12,-11,12,3,-2,1,8,3,12,-1,-10,10,12,7,6,2,4,12,10,-7,-4,2,7,3,11,8,9,-6,-5,-3,-9,12,6,-8,6,-2,-5,10,-8,-5,9,-9,1,9,-1,12,-6,7,10,2,-5,2,1,7,6,-8,-3,-3,8,-6,-5,3,8,2,12,0,9,-3,-1,12,5,-9,8,7,-7,-7,-12,3,12,-6,9,2,-10,-7,-10,11,-1,0,-12,-10,-2,3,-4,-3,-2,-4,6,-5,12,12,0,-3,-6,-8,-6,-6,-4,-8,5,10,10,10,1,-6,1,-8,10,3,12,-5,-8,8,8,-3,10,5,-4,3,-6,4,-10,12,-6,3,11,8,-6,-3,-1,-3,-8,12,3,11,7,12,-3,4,2,-8,-11,-11,11,1,-9,-6,-8,8,3,-1,11,12,3,0,4,-10,12,9,8,-10,12,10,12,0};
int BRIEF::gaussian_bit_pattern_31_y_b[256] = {5,-12,2,-13,12,6,-4,-8,-9,9,-9,12,6,0,-3,5,-1,12,-8,-8,1,-3,12,-2,-10,10,-3,7,11,-7,-1,-5,-13,12,4,7,-10,12,-13,2,3,-9,7,3,-10,0,1,12,-4,-12,-4,8,-7,-12,6,-10,5,12,8,7,8,-6,12,5,-13,5,-7,-11,-13,-1,2,12,6,-4,-3,12,5,4,2,1,5,-6,-7,-12,12,0,-13,9,-6,12,6,3,5,12,9,11,10,3,-6,-13,3,9,-6,-8,-4,-2,0,-8,3,-4,10,12,0,-6,-11,7,7,12,2,12,-8,-2,-13,0,-2,1,-4,-11,4,12,8,8,-13,12,7,-9,-8,9,-3,-12,0,12,-2,10,-4,-13,12,-6,3,-5,1,-11,-7,-5,6,6,1,-8,-8,9,3,7,-8,8,3,-9,-5,8,12,9,-5,11,-13,2,0,-10,-7,9,11,5,6,-2,7,-2,7,-13,-8,-9,5,10,-13,-13,-1,-9,-13,2,12,-10,-6,-6,-9,-7,-13,5,-13,-3,-12,-1,3,-9,1,-8,9,12,-5,7,-8,-12,5,9,5,4,3,12,11,-13,12,4,6,12,1,1,1,-13,-13,4,-2,-3,-2,10,-9,-1,-2,-8,5,10,5,5,11,-6,-12,9,4,-2,-2,-11};