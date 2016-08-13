

#ifndef BRIEF_H
#define BRIEF_H

#include "Types.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include "immintrin.h"

#ifndef N_DIM_BINARYDESCRIPTOR
#define N_DIM_BINARYDESCRIPTOR 256
#else
assert(256 == 256)
#endif

#ifndef SIZE_BITS_HAMING
#define SIZE_BITS_HAMING 64
#else
assert(64 == 64)
#endif

class BRIEF {
public:
    BRIEF(int n_rows_, int n_cols_);

template<class T>
void rbrief(unsigned char *image_src, const int height_image, const int width_image, const int n_channels,
             const int stride_image, T *AoS, const int n_features) {
    for (int j = 0; j < n_features; j++) {
        if ((AoS[j].x > diag_length_pattern) && AoS[j].x < (width_image - diag_length_pattern)
            && (AoS[j].y > diag_length_pattern) && AoS[j].y < (height_image - diag_length_pattern)) {
            float cos_angle = std::cos(AoS[j].angle);
            float sin_angle = std::sin(AoS[j].angle);
            unsigned char *image_center = image_src + AoS[j].y * stride_image + AoS[j].x * n_channels;
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
            f[0] = (*(image_center + i_y_a[0*32]*stride_image + i_x_a[0*32])
                        < *(image_center + i_y_b[0*32]*stride_image + i_x_b[0*32]));
            f[0] |= (*(image_center + i_y_a[0*32 + 1]*stride_image + i_x_a[0*32 + 1]) < *(image_center + i_y_b[0*32 + 1]*stride_image + i_x_b[0*32 + 1])) << 1;
            f[0] |= (*(image_center + i_y_a[0*32 + 2]*stride_image + i_x_a[0*32 + 2]) < *(image_center + i_y_b[0*32 + 2]*stride_image + i_x_b[0*32 + 2])) << 2;
            f[0] |= (*(image_center + i_y_a[0*32 + 3]*stride_image + i_x_a[0*32 + 3]) < *(image_center + i_y_b[0*32 + 3]*stride_image + i_x_b[0*32 + 3])) << 3;
            f[0] |= (*(image_center + i_y_a[0*32 + 4]*stride_image + i_x_a[0*32 + 4]) < *(image_center + i_y_b[0*32 + 4]*stride_image + i_x_b[0*32 + 4])) << 4;
            f[0] |= (*(image_center + i_y_a[0*32 + 5]*stride_image + i_x_a[0*32 + 5]) < *(image_center + i_y_b[0*32 + 5]*stride_image + i_x_b[0*32 + 5])) << 5;
            f[0] |= (*(image_center + i_y_a[0*32 + 6]*stride_image + i_x_a[0*32 + 6]) < *(image_center + i_y_b[0*32 + 6]*stride_image + i_x_b[0*32 + 6])) << 6;
            f[0] |= (*(image_center + i_y_a[0*32 + 7]*stride_image + i_x_a[0*32 + 7]) < *(image_center + i_y_b[0*32 + 7]*stride_image + i_x_b[0*32 + 7])) << 7;
            f[0] |= (*(image_center + i_y_a[0*32 + 8]*stride_image + i_x_a[0*32 + 8]) < *(image_center + i_y_b[0*32 + 8]*stride_image + i_x_b[0*32 + 8])) << 8;
            f[0] |= (*(image_center + i_y_a[0*32 + 9]*stride_image + i_x_a[0*32 + 9]) < *(image_center + i_y_b[0*32 + 9]*stride_image + i_x_b[0*32 + 9])) << 9;
            f[0] |= (*(image_center + i_y_a[0*32 + 10]*stride_image + i_x_a[0*32 + 10]) < *(image_center + i_y_b[0*32 + 10]*stride_image + i_x_b[0*32 + 10])) << 10;
            f[0] |= (*(image_center + i_y_a[0*32 + 11]*stride_image + i_x_a[0*32 + 11]) < *(image_center + i_y_b[0*32 + 11]*stride_image + i_x_b[0*32 + 11])) << 11;
            f[0] |= (*(image_center + i_y_a[0*32 + 12]*stride_image + i_x_a[0*32 + 12]) < *(image_center + i_y_b[0*32 + 12]*stride_image + i_x_b[0*32 + 12])) << 12;
            f[0] |= (*(image_center + i_y_a[0*32 + 13]*stride_image + i_x_a[0*32 + 13]) < *(image_center + i_y_b[0*32 + 13]*stride_image + i_x_b[0*32 + 13])) << 13;
            f[0] |= (*(image_center + i_y_a[0*32 + 14]*stride_image + i_x_a[0*32 + 14]) < *(image_center + i_y_b[0*32 + 14]*stride_image + i_x_b[0*32 + 14])) << 14;
            f[0] |= (*(image_center + i_y_a[0*32 + 15]*stride_image + i_x_a[0*32 + 15]) < *(image_center + i_y_b[0*32 + 15]*stride_image + i_x_b[0*32 + 15])) << 15;
            f[0] |= (*(image_center + i_y_a[0*32 + 16]*stride_image + i_x_a[0*32 + 16]) < *(image_center + i_y_b[0*32 + 16]*stride_image + i_x_b[0*32 + 16])) << 16;
            f[0] |= (*(image_center + i_y_a[0*32 + 17]*stride_image + i_x_a[0*32 + 17]) < *(image_center + i_y_b[0*32 + 17]*stride_image + i_x_b[0*32 + 17])) << 17;
            f[0] |= (*(image_center + i_y_a[0*32 + 18]*stride_image + i_x_a[0*32 + 18]) < *(image_center + i_y_b[0*32 + 18]*stride_image + i_x_b[0*32 + 18])) << 18;
            f[0] |= (*(image_center + i_y_a[0*32 + 19]*stride_image + i_x_a[0*32 + 19]) < *(image_center + i_y_b[0*32 + 19]*stride_image + i_x_b[0*32 + 19])) << 19;
            f[0] |= (*(image_center + i_y_a[0*32 + 20]*stride_image + i_x_a[0*32 + 20]) < *(image_center + i_y_b[0*32 + 20]*stride_image + i_x_b[0*32 + 20])) << 20;
            f[0] |= (*(image_center + i_y_a[0*32 + 21]*stride_image + i_x_a[0*32 + 21]) < *(image_center + i_y_b[0*32 + 21]*stride_image + i_x_b[0*32 + 21])) << 21;
            f[0] |= (*(image_center + i_y_a[0*32 + 22]*stride_image + i_x_a[0*32 + 22]) < *(image_center + i_y_b[0*32 + 22]*stride_image + i_x_b[0*32 + 22])) << 22;
            f[0] |= (*(image_center + i_y_a[0*32 + 23]*stride_image + i_x_a[0*32 + 23]) < *(image_center + i_y_b[0*32 + 23]*stride_image + i_x_b[0*32 + 23])) << 23;
            f[0] |= (*(image_center + i_y_a[0*32 + 24]*stride_image + i_x_a[0*32 + 24]) < *(image_center + i_y_b[0*32 + 24]*stride_image + i_x_b[0*32 + 24])) << 24;
            f[0] |= (*(image_center + i_y_a[0*32 + 25]*stride_image + i_x_a[0*32 + 25]) < *(image_center + i_y_b[0*32 + 25]*stride_image + i_x_b[0*32 + 25])) << 25;
            f[0] |= (*(image_center + i_y_a[0*32 + 26]*stride_image + i_x_a[0*32 + 26]) < *(image_center + i_y_b[0*32 + 26]*stride_image + i_x_b[0*32 + 26])) << 26;
            f[0] |= (*(image_center + i_y_a[0*32 + 27]*stride_image + i_x_a[0*32 + 27]) < *(image_center + i_y_b[0*32 + 27]*stride_image + i_x_b[0*32 + 27])) << 27;
            f[0] |= (*(image_center + i_y_a[0*32 + 28]*stride_image + i_x_a[0*32 + 28]) < *(image_center + i_y_b[0*32 + 28]*stride_image + i_x_b[0*32 + 28])) << 28;
            f[0] |= (*(image_center + i_y_a[0*32 + 29]*stride_image + i_x_a[0*32 + 29]) < *(image_center + i_y_b[0*32 + 29]*stride_image + i_x_b[0*32 + 29])) << 29;
            f[0] |= (*(image_center + i_y_a[0*32 + 30]*stride_image + i_x_a[0*32 + 30]) < *(image_center + i_y_b[0*32 + 30]*stride_image + i_x_b[0*32 + 30])) << 30;
            f[0] |= (*(image_center + i_y_a[0*32 + 31]*stride_image + i_x_a[0*32 + 31]) < *(image_center + i_y_b[0*32 + 31]*stride_image + i_x_b[0*32 + 31])) << 31;
            
            f[1] = (*(image_center + i_y_a[1*32]*stride_image + i_x_a[1*32])
                        < *(image_center + i_y_b[1*32]*stride_image + i_x_b[1*32]));
            f[1] |= (*(image_center + i_y_a[1*32 + 1]*stride_image + i_x_a[1*32 + 1]) < *(image_center + i_y_b[1*32 + 1]*stride_image + i_x_b[1*32 + 1])) << 1;
            f[1] |= (*(image_center + i_y_a[1*32 + 2]*stride_image + i_x_a[1*32 + 2]) < *(image_center + i_y_b[1*32 + 2]*stride_image + i_x_b[1*32 + 2])) << 2;
            f[1] |= (*(image_center + i_y_a[1*32 + 3]*stride_image + i_x_a[1*32 + 3]) < *(image_center + i_y_b[1*32 + 3]*stride_image + i_x_b[1*32 + 3])) << 3;
            f[1] |= (*(image_center + i_y_a[1*32 + 4]*stride_image + i_x_a[1*32 + 4]) < *(image_center + i_y_b[1*32 + 4]*stride_image + i_x_b[1*32 + 4])) << 4;
            f[1] |= (*(image_center + i_y_a[1*32 + 5]*stride_image + i_x_a[1*32 + 5]) < *(image_center + i_y_b[1*32 + 5]*stride_image + i_x_b[1*32 + 5])) << 5;
            f[1] |= (*(image_center + i_y_a[1*32 + 6]*stride_image + i_x_a[1*32 + 6]) < *(image_center + i_y_b[1*32 + 6]*stride_image + i_x_b[1*32 + 6])) << 6;
            f[1] |= (*(image_center + i_y_a[1*32 + 7]*stride_image + i_x_a[1*32 + 7]) < *(image_center + i_y_b[1*32 + 7]*stride_image + i_x_b[1*32 + 7])) << 7;
            f[1] |= (*(image_center + i_y_a[1*32 + 8]*stride_image + i_x_a[1*32 + 8]) < *(image_center + i_y_b[1*32 + 8]*stride_image + i_x_b[1*32 + 8])) << 8;
            f[1] |= (*(image_center + i_y_a[1*32 + 9]*stride_image + i_x_a[1*32 + 9]) < *(image_center + i_y_b[1*32 + 9]*stride_image + i_x_b[1*32 + 9])) << 9;
            f[1] |= (*(image_center + i_y_a[1*32 + 10]*stride_image + i_x_a[1*32 + 10]) < *(image_center + i_y_b[1*32 + 10]*stride_image + i_x_b[1*32 + 10])) << 10;
            f[1] |= (*(image_center + i_y_a[1*32 + 11]*stride_image + i_x_a[1*32 + 11]) < *(image_center + i_y_b[1*32 + 11]*stride_image + i_x_b[1*32 + 11])) << 11;
            f[1] |= (*(image_center + i_y_a[1*32 + 12]*stride_image + i_x_a[1*32 + 12]) < *(image_center + i_y_b[1*32 + 12]*stride_image + i_x_b[1*32 + 12])) << 12;
            f[1] |= (*(image_center + i_y_a[1*32 + 13]*stride_image + i_x_a[1*32 + 13]) < *(image_center + i_y_b[1*32 + 13]*stride_image + i_x_b[1*32 + 13])) << 13;
            f[1] |= (*(image_center + i_y_a[1*32 + 14]*stride_image + i_x_a[1*32 + 14]) < *(image_center + i_y_b[1*32 + 14]*stride_image + i_x_b[1*32 + 14])) << 14;
            f[1] |= (*(image_center + i_y_a[1*32 + 15]*stride_image + i_x_a[1*32 + 15]) < *(image_center + i_y_b[1*32 + 15]*stride_image + i_x_b[1*32 + 15])) << 15;
            f[1] |= (*(image_center + i_y_a[1*32 + 16]*stride_image + i_x_a[1*32 + 16]) < *(image_center + i_y_b[1*32 + 16]*stride_image + i_x_b[1*32 + 16])) << 16;
            f[1] |= (*(image_center + i_y_a[1*32 + 17]*stride_image + i_x_a[1*32 + 17]) < *(image_center + i_y_b[1*32 + 17]*stride_image + i_x_b[1*32 + 17])) << 17;
            f[1] |= (*(image_center + i_y_a[1*32 + 18]*stride_image + i_x_a[1*32 + 18]) < *(image_center + i_y_b[1*32 + 18]*stride_image + i_x_b[1*32 + 18])) << 18;
            f[1] |= (*(image_center + i_y_a[1*32 + 19]*stride_image + i_x_a[1*32 + 19]) < *(image_center + i_y_b[1*32 + 19]*stride_image + i_x_b[1*32 + 19])) << 19;
            f[1] |= (*(image_center + i_y_a[1*32 + 20]*stride_image + i_x_a[1*32 + 20]) < *(image_center + i_y_b[1*32 + 20]*stride_image + i_x_b[1*32 + 20])) << 20;
            f[1] |= (*(image_center + i_y_a[1*32 + 21]*stride_image + i_x_a[1*32 + 21]) < *(image_center + i_y_b[1*32 + 21]*stride_image + i_x_b[1*32 + 21])) << 21;
            f[1] |= (*(image_center + i_y_a[1*32 + 22]*stride_image + i_x_a[1*32 + 22]) < *(image_center + i_y_b[1*32 + 22]*stride_image + i_x_b[1*32 + 22])) << 22;
            f[1] |= (*(image_center + i_y_a[1*32 + 23]*stride_image + i_x_a[1*32 + 23]) < *(image_center + i_y_b[1*32 + 23]*stride_image + i_x_b[1*32 + 23])) << 23;
            f[1] |= (*(image_center + i_y_a[1*32 + 24]*stride_image + i_x_a[1*32 + 24]) < *(image_center + i_y_b[1*32 + 24]*stride_image + i_x_b[1*32 + 24])) << 24;
            f[1] |= (*(image_center + i_y_a[1*32 + 25]*stride_image + i_x_a[1*32 + 25]) < *(image_center + i_y_b[1*32 + 25]*stride_image + i_x_b[1*32 + 25])) << 25;
            f[1] |= (*(image_center + i_y_a[1*32 + 26]*stride_image + i_x_a[1*32 + 26]) < *(image_center + i_y_b[1*32 + 26]*stride_image + i_x_b[1*32 + 26])) << 26;
            f[1] |= (*(image_center + i_y_a[1*32 + 27]*stride_image + i_x_a[1*32 + 27]) < *(image_center + i_y_b[1*32 + 27]*stride_image + i_x_b[1*32 + 27])) << 27;
            f[1] |= (*(image_center + i_y_a[1*32 + 28]*stride_image + i_x_a[1*32 + 28]) < *(image_center + i_y_b[1*32 + 28]*stride_image + i_x_b[1*32 + 28])) << 28;
            f[1] |= (*(image_center + i_y_a[1*32 + 29]*stride_image + i_x_a[1*32 + 29]) < *(image_center + i_y_b[1*32 + 29]*stride_image + i_x_b[1*32 + 29])) << 29;
            f[1] |= (*(image_center + i_y_a[1*32 + 30]*stride_image + i_x_a[1*32 + 30]) < *(image_center + i_y_b[1*32 + 30]*stride_image + i_x_b[1*32 + 30])) << 30;
            f[1] |= (*(image_center + i_y_a[1*32 + 31]*stride_image + i_x_a[1*32 + 31]) < *(image_center + i_y_b[1*32 + 31]*stride_image + i_x_b[1*32 + 31])) << 31;
            
            f[2] = (*(image_center + i_y_a[2*32]*stride_image + i_x_a[2*32])
                        < *(image_center + i_y_b[2*32]*stride_image + i_x_b[2*32]));
            f[2] |= (*(image_center + i_y_a[2*32 + 1]*stride_image + i_x_a[2*32 + 1]) < *(image_center + i_y_b[2*32 + 1]*stride_image + i_x_b[2*32 + 1])) << 1;
            f[2] |= (*(image_center + i_y_a[2*32 + 2]*stride_image + i_x_a[2*32 + 2]) < *(image_center + i_y_b[2*32 + 2]*stride_image + i_x_b[2*32 + 2])) << 2;
            f[2] |= (*(image_center + i_y_a[2*32 + 3]*stride_image + i_x_a[2*32 + 3]) < *(image_center + i_y_b[2*32 + 3]*stride_image + i_x_b[2*32 + 3])) << 3;
            f[2] |= (*(image_center + i_y_a[2*32 + 4]*stride_image + i_x_a[2*32 + 4]) < *(image_center + i_y_b[2*32 + 4]*stride_image + i_x_b[2*32 + 4])) << 4;
            f[2] |= (*(image_center + i_y_a[2*32 + 5]*stride_image + i_x_a[2*32 + 5]) < *(image_center + i_y_b[2*32 + 5]*stride_image + i_x_b[2*32 + 5])) << 5;
            f[2] |= (*(image_center + i_y_a[2*32 + 6]*stride_image + i_x_a[2*32 + 6]) < *(image_center + i_y_b[2*32 + 6]*stride_image + i_x_b[2*32 + 6])) << 6;
            f[2] |= (*(image_center + i_y_a[2*32 + 7]*stride_image + i_x_a[2*32 + 7]) < *(image_center + i_y_b[2*32 + 7]*stride_image + i_x_b[2*32 + 7])) << 7;
            f[2] |= (*(image_center + i_y_a[2*32 + 8]*stride_image + i_x_a[2*32 + 8]) < *(image_center + i_y_b[2*32 + 8]*stride_image + i_x_b[2*32 + 8])) << 8;
            f[2] |= (*(image_center + i_y_a[2*32 + 9]*stride_image + i_x_a[2*32 + 9]) < *(image_center + i_y_b[2*32 + 9]*stride_image + i_x_b[2*32 + 9])) << 9;
            f[2] |= (*(image_center + i_y_a[2*32 + 10]*stride_image + i_x_a[2*32 + 10]) < *(image_center + i_y_b[2*32 + 10]*stride_image + i_x_b[2*32 + 10])) << 10;
            f[2] |= (*(image_center + i_y_a[2*32 + 11]*stride_image + i_x_a[2*32 + 11]) < *(image_center + i_y_b[2*32 + 11]*stride_image + i_x_b[2*32 + 11])) << 11;
            f[2] |= (*(image_center + i_y_a[2*32 + 12]*stride_image + i_x_a[2*32 + 12]) < *(image_center + i_y_b[2*32 + 12]*stride_image + i_x_b[2*32 + 12])) << 12;
            f[2] |= (*(image_center + i_y_a[2*32 + 13]*stride_image + i_x_a[2*32 + 13]) < *(image_center + i_y_b[2*32 + 13]*stride_image + i_x_b[2*32 + 13])) << 13;
            f[2] |= (*(image_center + i_y_a[2*32 + 14]*stride_image + i_x_a[2*32 + 14]) < *(image_center + i_y_b[2*32 + 14]*stride_image + i_x_b[2*32 + 14])) << 14;
            f[2] |= (*(image_center + i_y_a[2*32 + 15]*stride_image + i_x_a[2*32 + 15]) < *(image_center + i_y_b[2*32 + 15]*stride_image + i_x_b[2*32 + 15])) << 15;
            f[2] |= (*(image_center + i_y_a[2*32 + 16]*stride_image + i_x_a[2*32 + 16]) < *(image_center + i_y_b[2*32 + 16]*stride_image + i_x_b[2*32 + 16])) << 16;
            f[2] |= (*(image_center + i_y_a[2*32 + 17]*stride_image + i_x_a[2*32 + 17]) < *(image_center + i_y_b[2*32 + 17]*stride_image + i_x_b[2*32 + 17])) << 17;
            f[2] |= (*(image_center + i_y_a[2*32 + 18]*stride_image + i_x_a[2*32 + 18]) < *(image_center + i_y_b[2*32 + 18]*stride_image + i_x_b[2*32 + 18])) << 18;
            f[2] |= (*(image_center + i_y_a[2*32 + 19]*stride_image + i_x_a[2*32 + 19]) < *(image_center + i_y_b[2*32 + 19]*stride_image + i_x_b[2*32 + 19])) << 19;
            f[2] |= (*(image_center + i_y_a[2*32 + 20]*stride_image + i_x_a[2*32 + 20]) < *(image_center + i_y_b[2*32 + 20]*stride_image + i_x_b[2*32 + 20])) << 20;
            f[2] |= (*(image_center + i_y_a[2*32 + 21]*stride_image + i_x_a[2*32 + 21]) < *(image_center + i_y_b[2*32 + 21]*stride_image + i_x_b[2*32 + 21])) << 21;
            f[2] |= (*(image_center + i_y_a[2*32 + 22]*stride_image + i_x_a[2*32 + 22]) < *(image_center + i_y_b[2*32 + 22]*stride_image + i_x_b[2*32 + 22])) << 22;
            f[2] |= (*(image_center + i_y_a[2*32 + 23]*stride_image + i_x_a[2*32 + 23]) < *(image_center + i_y_b[2*32 + 23]*stride_image + i_x_b[2*32 + 23])) << 23;
            f[2] |= (*(image_center + i_y_a[2*32 + 24]*stride_image + i_x_a[2*32 + 24]) < *(image_center + i_y_b[2*32 + 24]*stride_image + i_x_b[2*32 + 24])) << 24;
            f[2] |= (*(image_center + i_y_a[2*32 + 25]*stride_image + i_x_a[2*32 + 25]) < *(image_center + i_y_b[2*32 + 25]*stride_image + i_x_b[2*32 + 25])) << 25;
            f[2] |= (*(image_center + i_y_a[2*32 + 26]*stride_image + i_x_a[2*32 + 26]) < *(image_center + i_y_b[2*32 + 26]*stride_image + i_x_b[2*32 + 26])) << 26;
            f[2] |= (*(image_center + i_y_a[2*32 + 27]*stride_image + i_x_a[2*32 + 27]) < *(image_center + i_y_b[2*32 + 27]*stride_image + i_x_b[2*32 + 27])) << 27;
            f[2] |= (*(image_center + i_y_a[2*32 + 28]*stride_image + i_x_a[2*32 + 28]) < *(image_center + i_y_b[2*32 + 28]*stride_image + i_x_b[2*32 + 28])) << 28;
            f[2] |= (*(image_center + i_y_a[2*32 + 29]*stride_image + i_x_a[2*32 + 29]) < *(image_center + i_y_b[2*32 + 29]*stride_image + i_x_b[2*32 + 29])) << 29;
            f[2] |= (*(image_center + i_y_a[2*32 + 30]*stride_image + i_x_a[2*32 + 30]) < *(image_center + i_y_b[2*32 + 30]*stride_image + i_x_b[2*32 + 30])) << 30;
            f[2] |= (*(image_center + i_y_a[2*32 + 31]*stride_image + i_x_a[2*32 + 31]) < *(image_center + i_y_b[2*32 + 31]*stride_image + i_x_b[2*32 + 31])) << 31;
            
            f[3] = (*(image_center + i_y_a[3*32]*stride_image + i_x_a[3*32])
                        < *(image_center + i_y_b[3*32]*stride_image + i_x_b[3*32]));
            f[3] |= (*(image_center + i_y_a[3*32 + 1]*stride_image + i_x_a[3*32 + 1]) < *(image_center + i_y_b[3*32 + 1]*stride_image + i_x_b[3*32 + 1])) << 1;
            f[3] |= (*(image_center + i_y_a[3*32 + 2]*stride_image + i_x_a[3*32 + 2]) < *(image_center + i_y_b[3*32 + 2]*stride_image + i_x_b[3*32 + 2])) << 2;
            f[3] |= (*(image_center + i_y_a[3*32 + 3]*stride_image + i_x_a[3*32 + 3]) < *(image_center + i_y_b[3*32 + 3]*stride_image + i_x_b[3*32 + 3])) << 3;
            f[3] |= (*(image_center + i_y_a[3*32 + 4]*stride_image + i_x_a[3*32 + 4]) < *(image_center + i_y_b[3*32 + 4]*stride_image + i_x_b[3*32 + 4])) << 4;
            f[3] |= (*(image_center + i_y_a[3*32 + 5]*stride_image + i_x_a[3*32 + 5]) < *(image_center + i_y_b[3*32 + 5]*stride_image + i_x_b[3*32 + 5])) << 5;
            f[3] |= (*(image_center + i_y_a[3*32 + 6]*stride_image + i_x_a[3*32 + 6]) < *(image_center + i_y_b[3*32 + 6]*stride_image + i_x_b[3*32 + 6])) << 6;
            f[3] |= (*(image_center + i_y_a[3*32 + 7]*stride_image + i_x_a[3*32 + 7]) < *(image_center + i_y_b[3*32 + 7]*stride_image + i_x_b[3*32 + 7])) << 7;
            f[3] |= (*(image_center + i_y_a[3*32 + 8]*stride_image + i_x_a[3*32 + 8]) < *(image_center + i_y_b[3*32 + 8]*stride_image + i_x_b[3*32 + 8])) << 8;
            f[3] |= (*(image_center + i_y_a[3*32 + 9]*stride_image + i_x_a[3*32 + 9]) < *(image_center + i_y_b[3*32 + 9]*stride_image + i_x_b[3*32 + 9])) << 9;
            f[3] |= (*(image_center + i_y_a[3*32 + 10]*stride_image + i_x_a[3*32 + 10]) < *(image_center + i_y_b[3*32 + 10]*stride_image + i_x_b[3*32 + 10])) << 10;
            f[3] |= (*(image_center + i_y_a[3*32 + 11]*stride_image + i_x_a[3*32 + 11]) < *(image_center + i_y_b[3*32 + 11]*stride_image + i_x_b[3*32 + 11])) << 11;
            f[3] |= (*(image_center + i_y_a[3*32 + 12]*stride_image + i_x_a[3*32 + 12]) < *(image_center + i_y_b[3*32 + 12]*stride_image + i_x_b[3*32 + 12])) << 12;
            f[3] |= (*(image_center + i_y_a[3*32 + 13]*stride_image + i_x_a[3*32 + 13]) < *(image_center + i_y_b[3*32 + 13]*stride_image + i_x_b[3*32 + 13])) << 13;
            f[3] |= (*(image_center + i_y_a[3*32 + 14]*stride_image + i_x_a[3*32 + 14]) < *(image_center + i_y_b[3*32 + 14]*stride_image + i_x_b[3*32 + 14])) << 14;
            f[3] |= (*(image_center + i_y_a[3*32 + 15]*stride_image + i_x_a[3*32 + 15]) < *(image_center + i_y_b[3*32 + 15]*stride_image + i_x_b[3*32 + 15])) << 15;
            f[3] |= (*(image_center + i_y_a[3*32 + 16]*stride_image + i_x_a[3*32 + 16]) < *(image_center + i_y_b[3*32 + 16]*stride_image + i_x_b[3*32 + 16])) << 16;
            f[3] |= (*(image_center + i_y_a[3*32 + 17]*stride_image + i_x_a[3*32 + 17]) < *(image_center + i_y_b[3*32 + 17]*stride_image + i_x_b[3*32 + 17])) << 17;
            f[3] |= (*(image_center + i_y_a[3*32 + 18]*stride_image + i_x_a[3*32 + 18]) < *(image_center + i_y_b[3*32 + 18]*stride_image + i_x_b[3*32 + 18])) << 18;
            f[3] |= (*(image_center + i_y_a[3*32 + 19]*stride_image + i_x_a[3*32 + 19]) < *(image_center + i_y_b[3*32 + 19]*stride_image + i_x_b[3*32 + 19])) << 19;
            f[3] |= (*(image_center + i_y_a[3*32 + 20]*stride_image + i_x_a[3*32 + 20]) < *(image_center + i_y_b[3*32 + 20]*stride_image + i_x_b[3*32 + 20])) << 20;
            f[3] |= (*(image_center + i_y_a[3*32 + 21]*stride_image + i_x_a[3*32 + 21]) < *(image_center + i_y_b[3*32 + 21]*stride_image + i_x_b[3*32 + 21])) << 21;
            f[3] |= (*(image_center + i_y_a[3*32 + 22]*stride_image + i_x_a[3*32 + 22]) < *(image_center + i_y_b[3*32 + 22]*stride_image + i_x_b[3*32 + 22])) << 22;
            f[3] |= (*(image_center + i_y_a[3*32 + 23]*stride_image + i_x_a[3*32 + 23]) < *(image_center + i_y_b[3*32 + 23]*stride_image + i_x_b[3*32 + 23])) << 23;
            f[3] |= (*(image_center + i_y_a[3*32 + 24]*stride_image + i_x_a[3*32 + 24]) < *(image_center + i_y_b[3*32 + 24]*stride_image + i_x_b[3*32 + 24])) << 24;
            f[3] |= (*(image_center + i_y_a[3*32 + 25]*stride_image + i_x_a[3*32 + 25]) < *(image_center + i_y_b[3*32 + 25]*stride_image + i_x_b[3*32 + 25])) << 25;
            f[3] |= (*(image_center + i_y_a[3*32 + 26]*stride_image + i_x_a[3*32 + 26]) < *(image_center + i_y_b[3*32 + 26]*stride_image + i_x_b[3*32 + 26])) << 26;
            f[3] |= (*(image_center + i_y_a[3*32 + 27]*stride_image + i_x_a[3*32 + 27]) < *(image_center + i_y_b[3*32 + 27]*stride_image + i_x_b[3*32 + 27])) << 27;
            f[3] |= (*(image_center + i_y_a[3*32 + 28]*stride_image + i_x_a[3*32 + 28]) < *(image_center + i_y_b[3*32 + 28]*stride_image + i_x_b[3*32 + 28])) << 28;
            f[3] |= (*(image_center + i_y_a[3*32 + 29]*stride_image + i_x_a[3*32 + 29]) < *(image_center + i_y_b[3*32 + 29]*stride_image + i_x_b[3*32 + 29])) << 29;
            f[3] |= (*(image_center + i_y_a[3*32 + 30]*stride_image + i_x_a[3*32 + 30]) < *(image_center + i_y_b[3*32 + 30]*stride_image + i_x_b[3*32 + 30])) << 30;
            f[3] |= (*(image_center + i_y_a[3*32 + 31]*stride_image + i_x_a[3*32 + 31]) < *(image_center + i_y_b[3*32 + 31]*stride_image + i_x_b[3*32 + 31])) << 31;
            
            f[4] = (*(image_center + i_y_a[4*32]*stride_image + i_x_a[4*32])
                        < *(image_center + i_y_b[4*32]*stride_image + i_x_b[4*32]));
            f[4] |= (*(image_center + i_y_a[4*32 + 1]*stride_image + i_x_a[4*32 + 1]) < *(image_center + i_y_b[4*32 + 1]*stride_image + i_x_b[4*32 + 1])) << 1;
            f[4] |= (*(image_center + i_y_a[4*32 + 2]*stride_image + i_x_a[4*32 + 2]) < *(image_center + i_y_b[4*32 + 2]*stride_image + i_x_b[4*32 + 2])) << 2;
            f[4] |= (*(image_center + i_y_a[4*32 + 3]*stride_image + i_x_a[4*32 + 3]) < *(image_center + i_y_b[4*32 + 3]*stride_image + i_x_b[4*32 + 3])) << 3;
            f[4] |= (*(image_center + i_y_a[4*32 + 4]*stride_image + i_x_a[4*32 + 4]) < *(image_center + i_y_b[4*32 + 4]*stride_image + i_x_b[4*32 + 4])) << 4;
            f[4] |= (*(image_center + i_y_a[4*32 + 5]*stride_image + i_x_a[4*32 + 5]) < *(image_center + i_y_b[4*32 + 5]*stride_image + i_x_b[4*32 + 5])) << 5;
            f[4] |= (*(image_center + i_y_a[4*32 + 6]*stride_image + i_x_a[4*32 + 6]) < *(image_center + i_y_b[4*32 + 6]*stride_image + i_x_b[4*32 + 6])) << 6;
            f[4] |= (*(image_center + i_y_a[4*32 + 7]*stride_image + i_x_a[4*32 + 7]) < *(image_center + i_y_b[4*32 + 7]*stride_image + i_x_b[4*32 + 7])) << 7;
            f[4] |= (*(image_center + i_y_a[4*32 + 8]*stride_image + i_x_a[4*32 + 8]) < *(image_center + i_y_b[4*32 + 8]*stride_image + i_x_b[4*32 + 8])) << 8;
            f[4] |= (*(image_center + i_y_a[4*32 + 9]*stride_image + i_x_a[4*32 + 9]) < *(image_center + i_y_b[4*32 + 9]*stride_image + i_x_b[4*32 + 9])) << 9;
            f[4] |= (*(image_center + i_y_a[4*32 + 10]*stride_image + i_x_a[4*32 + 10]) < *(image_center + i_y_b[4*32 + 10]*stride_image + i_x_b[4*32 + 10])) << 10;
            f[4] |= (*(image_center + i_y_a[4*32 + 11]*stride_image + i_x_a[4*32 + 11]) < *(image_center + i_y_b[4*32 + 11]*stride_image + i_x_b[4*32 + 11])) << 11;
            f[4] |= (*(image_center + i_y_a[4*32 + 12]*stride_image + i_x_a[4*32 + 12]) < *(image_center + i_y_b[4*32 + 12]*stride_image + i_x_b[4*32 + 12])) << 12;
            f[4] |= (*(image_center + i_y_a[4*32 + 13]*stride_image + i_x_a[4*32 + 13]) < *(image_center + i_y_b[4*32 + 13]*stride_image + i_x_b[4*32 + 13])) << 13;
            f[4] |= (*(image_center + i_y_a[4*32 + 14]*stride_image + i_x_a[4*32 + 14]) < *(image_center + i_y_b[4*32 + 14]*stride_image + i_x_b[4*32 + 14])) << 14;
            f[4] |= (*(image_center + i_y_a[4*32 + 15]*stride_image + i_x_a[4*32 + 15]) < *(image_center + i_y_b[4*32 + 15]*stride_image + i_x_b[4*32 + 15])) << 15;
            f[4] |= (*(image_center + i_y_a[4*32 + 16]*stride_image + i_x_a[4*32 + 16]) < *(image_center + i_y_b[4*32 + 16]*stride_image + i_x_b[4*32 + 16])) << 16;
            f[4] |= (*(image_center + i_y_a[4*32 + 17]*stride_image + i_x_a[4*32 + 17]) < *(image_center + i_y_b[4*32 + 17]*stride_image + i_x_b[4*32 + 17])) << 17;
            f[4] |= (*(image_center + i_y_a[4*32 + 18]*stride_image + i_x_a[4*32 + 18]) < *(image_center + i_y_b[4*32 + 18]*stride_image + i_x_b[4*32 + 18])) << 18;
            f[4] |= (*(image_center + i_y_a[4*32 + 19]*stride_image + i_x_a[4*32 + 19]) < *(image_center + i_y_b[4*32 + 19]*stride_image + i_x_b[4*32 + 19])) << 19;
            f[4] |= (*(image_center + i_y_a[4*32 + 20]*stride_image + i_x_a[4*32 + 20]) < *(image_center + i_y_b[4*32 + 20]*stride_image + i_x_b[4*32 + 20])) << 20;
            f[4] |= (*(image_center + i_y_a[4*32 + 21]*stride_image + i_x_a[4*32 + 21]) < *(image_center + i_y_b[4*32 + 21]*stride_image + i_x_b[4*32 + 21])) << 21;
            f[4] |= (*(image_center + i_y_a[4*32 + 22]*stride_image + i_x_a[4*32 + 22]) < *(image_center + i_y_b[4*32 + 22]*stride_image + i_x_b[4*32 + 22])) << 22;
            f[4] |= (*(image_center + i_y_a[4*32 + 23]*stride_image + i_x_a[4*32 + 23]) < *(image_center + i_y_b[4*32 + 23]*stride_image + i_x_b[4*32 + 23])) << 23;
            f[4] |= (*(image_center + i_y_a[4*32 + 24]*stride_image + i_x_a[4*32 + 24]) < *(image_center + i_y_b[4*32 + 24]*stride_image + i_x_b[4*32 + 24])) << 24;
            f[4] |= (*(image_center + i_y_a[4*32 + 25]*stride_image + i_x_a[4*32 + 25]) < *(image_center + i_y_b[4*32 + 25]*stride_image + i_x_b[4*32 + 25])) << 25;
            f[4] |= (*(image_center + i_y_a[4*32 + 26]*stride_image + i_x_a[4*32 + 26]) < *(image_center + i_y_b[4*32 + 26]*stride_image + i_x_b[4*32 + 26])) << 26;
            f[4] |= (*(image_center + i_y_a[4*32 + 27]*stride_image + i_x_a[4*32 + 27]) < *(image_center + i_y_b[4*32 + 27]*stride_image + i_x_b[4*32 + 27])) << 27;
            f[4] |= (*(image_center + i_y_a[4*32 + 28]*stride_image + i_x_a[4*32 + 28]) < *(image_center + i_y_b[4*32 + 28]*stride_image + i_x_b[4*32 + 28])) << 28;
            f[4] |= (*(image_center + i_y_a[4*32 + 29]*stride_image + i_x_a[4*32 + 29]) < *(image_center + i_y_b[4*32 + 29]*stride_image + i_x_b[4*32 + 29])) << 29;
            f[4] |= (*(image_center + i_y_a[4*32 + 30]*stride_image + i_x_a[4*32 + 30]) < *(image_center + i_y_b[4*32 + 30]*stride_image + i_x_b[4*32 + 30])) << 30;
            f[4] |= (*(image_center + i_y_a[4*32 + 31]*stride_image + i_x_a[4*32 + 31]) < *(image_center + i_y_b[4*32 + 31]*stride_image + i_x_b[4*32 + 31])) << 31;
            
            f[5] = (*(image_center + i_y_a[5*32]*stride_image + i_x_a[5*32])
                        < *(image_center + i_y_b[5*32]*stride_image + i_x_b[5*32]));
            f[5] |= (*(image_center + i_y_a[5*32 + 1]*stride_image + i_x_a[5*32 + 1]) < *(image_center + i_y_b[5*32 + 1]*stride_image + i_x_b[5*32 + 1])) << 1;
            f[5] |= (*(image_center + i_y_a[5*32 + 2]*stride_image + i_x_a[5*32 + 2]) < *(image_center + i_y_b[5*32 + 2]*stride_image + i_x_b[5*32 + 2])) << 2;
            f[5] |= (*(image_center + i_y_a[5*32 + 3]*stride_image + i_x_a[5*32 + 3]) < *(image_center + i_y_b[5*32 + 3]*stride_image + i_x_b[5*32 + 3])) << 3;
            f[5] |= (*(image_center + i_y_a[5*32 + 4]*stride_image + i_x_a[5*32 + 4]) < *(image_center + i_y_b[5*32 + 4]*stride_image + i_x_b[5*32 + 4])) << 4;
            f[5] |= (*(image_center + i_y_a[5*32 + 5]*stride_image + i_x_a[5*32 + 5]) < *(image_center + i_y_b[5*32 + 5]*stride_image + i_x_b[5*32 + 5])) << 5;
            f[5] |= (*(image_center + i_y_a[5*32 + 6]*stride_image + i_x_a[5*32 + 6]) < *(image_center + i_y_b[5*32 + 6]*stride_image + i_x_b[5*32 + 6])) << 6;
            f[5] |= (*(image_center + i_y_a[5*32 + 7]*stride_image + i_x_a[5*32 + 7]) < *(image_center + i_y_b[5*32 + 7]*stride_image + i_x_b[5*32 + 7])) << 7;
            f[5] |= (*(image_center + i_y_a[5*32 + 8]*stride_image + i_x_a[5*32 + 8]) < *(image_center + i_y_b[5*32 + 8]*stride_image + i_x_b[5*32 + 8])) << 8;
            f[5] |= (*(image_center + i_y_a[5*32 + 9]*stride_image + i_x_a[5*32 + 9]) < *(image_center + i_y_b[5*32 + 9]*stride_image + i_x_b[5*32 + 9])) << 9;
            f[5] |= (*(image_center + i_y_a[5*32 + 10]*stride_image + i_x_a[5*32 + 10]) < *(image_center + i_y_b[5*32 + 10]*stride_image + i_x_b[5*32 + 10])) << 10;
            f[5] |= (*(image_center + i_y_a[5*32 + 11]*stride_image + i_x_a[5*32 + 11]) < *(image_center + i_y_b[5*32 + 11]*stride_image + i_x_b[5*32 + 11])) << 11;
            f[5] |= (*(image_center + i_y_a[5*32 + 12]*stride_image + i_x_a[5*32 + 12]) < *(image_center + i_y_b[5*32 + 12]*stride_image + i_x_b[5*32 + 12])) << 12;
            f[5] |= (*(image_center + i_y_a[5*32 + 13]*stride_image + i_x_a[5*32 + 13]) < *(image_center + i_y_b[5*32 + 13]*stride_image + i_x_b[5*32 + 13])) << 13;
            f[5] |= (*(image_center + i_y_a[5*32 + 14]*stride_image + i_x_a[5*32 + 14]) < *(image_center + i_y_b[5*32 + 14]*stride_image + i_x_b[5*32 + 14])) << 14;
            f[5] |= (*(image_center + i_y_a[5*32 + 15]*stride_image + i_x_a[5*32 + 15]) < *(image_center + i_y_b[5*32 + 15]*stride_image + i_x_b[5*32 + 15])) << 15;
            f[5] |= (*(image_center + i_y_a[5*32 + 16]*stride_image + i_x_a[5*32 + 16]) < *(image_center + i_y_b[5*32 + 16]*stride_image + i_x_b[5*32 + 16])) << 16;
            f[5] |= (*(image_center + i_y_a[5*32 + 17]*stride_image + i_x_a[5*32 + 17]) < *(image_center + i_y_b[5*32 + 17]*stride_image + i_x_b[5*32 + 17])) << 17;
            f[5] |= (*(image_center + i_y_a[5*32 + 18]*stride_image + i_x_a[5*32 + 18]) < *(image_center + i_y_b[5*32 + 18]*stride_image + i_x_b[5*32 + 18])) << 18;
            f[5] |= (*(image_center + i_y_a[5*32 + 19]*stride_image + i_x_a[5*32 + 19]) < *(image_center + i_y_b[5*32 + 19]*stride_image + i_x_b[5*32 + 19])) << 19;
            f[5] |= (*(image_center + i_y_a[5*32 + 20]*stride_image + i_x_a[5*32 + 20]) < *(image_center + i_y_b[5*32 + 20]*stride_image + i_x_b[5*32 + 20])) << 20;
            f[5] |= (*(image_center + i_y_a[5*32 + 21]*stride_image + i_x_a[5*32 + 21]) < *(image_center + i_y_b[5*32 + 21]*stride_image + i_x_b[5*32 + 21])) << 21;
            f[5] |= (*(image_center + i_y_a[5*32 + 22]*stride_image + i_x_a[5*32 + 22]) < *(image_center + i_y_b[5*32 + 22]*stride_image + i_x_b[5*32 + 22])) << 22;
            f[5] |= (*(image_center + i_y_a[5*32 + 23]*stride_image + i_x_a[5*32 + 23]) < *(image_center + i_y_b[5*32 + 23]*stride_image + i_x_b[5*32 + 23])) << 23;
            f[5] |= (*(image_center + i_y_a[5*32 + 24]*stride_image + i_x_a[5*32 + 24]) < *(image_center + i_y_b[5*32 + 24]*stride_image + i_x_b[5*32 + 24])) << 24;
            f[5] |= (*(image_center + i_y_a[5*32 + 25]*stride_image + i_x_a[5*32 + 25]) < *(image_center + i_y_b[5*32 + 25]*stride_image + i_x_b[5*32 + 25])) << 25;
            f[5] |= (*(image_center + i_y_a[5*32 + 26]*stride_image + i_x_a[5*32 + 26]) < *(image_center + i_y_b[5*32 + 26]*stride_image + i_x_b[5*32 + 26])) << 26;
            f[5] |= (*(image_center + i_y_a[5*32 + 27]*stride_image + i_x_a[5*32 + 27]) < *(image_center + i_y_b[5*32 + 27]*stride_image + i_x_b[5*32 + 27])) << 27;
            f[5] |= (*(image_center + i_y_a[5*32 + 28]*stride_image + i_x_a[5*32 + 28]) < *(image_center + i_y_b[5*32 + 28]*stride_image + i_x_b[5*32 + 28])) << 28;
            f[5] |= (*(image_center + i_y_a[5*32 + 29]*stride_image + i_x_a[5*32 + 29]) < *(image_center + i_y_b[5*32 + 29]*stride_image + i_x_b[5*32 + 29])) << 29;
            f[5] |= (*(image_center + i_y_a[5*32 + 30]*stride_image + i_x_a[5*32 + 30]) < *(image_center + i_y_b[5*32 + 30]*stride_image + i_x_b[5*32 + 30])) << 30;
            f[5] |= (*(image_center + i_y_a[5*32 + 31]*stride_image + i_x_a[5*32 + 31]) < *(image_center + i_y_b[5*32 + 31]*stride_image + i_x_b[5*32 + 31])) << 31;
            
            f[6] = (*(image_center + i_y_a[6*32]*stride_image + i_x_a[6*32])
                        < *(image_center + i_y_b[6*32]*stride_image + i_x_b[6*32]));
            f[6] |= (*(image_center + i_y_a[6*32 + 1]*stride_image + i_x_a[6*32 + 1]) < *(image_center + i_y_b[6*32 + 1]*stride_image + i_x_b[6*32 + 1])) << 1;
            f[6] |= (*(image_center + i_y_a[6*32 + 2]*stride_image + i_x_a[6*32 + 2]) < *(image_center + i_y_b[6*32 + 2]*stride_image + i_x_b[6*32 + 2])) << 2;
            f[6] |= (*(image_center + i_y_a[6*32 + 3]*stride_image + i_x_a[6*32 + 3]) < *(image_center + i_y_b[6*32 + 3]*stride_image + i_x_b[6*32 + 3])) << 3;
            f[6] |= (*(image_center + i_y_a[6*32 + 4]*stride_image + i_x_a[6*32 + 4]) < *(image_center + i_y_b[6*32 + 4]*stride_image + i_x_b[6*32 + 4])) << 4;
            f[6] |= (*(image_center + i_y_a[6*32 + 5]*stride_image + i_x_a[6*32 + 5]) < *(image_center + i_y_b[6*32 + 5]*stride_image + i_x_b[6*32 + 5])) << 5;
            f[6] |= (*(image_center + i_y_a[6*32 + 6]*stride_image + i_x_a[6*32 + 6]) < *(image_center + i_y_b[6*32 + 6]*stride_image + i_x_b[6*32 + 6])) << 6;
            f[6] |= (*(image_center + i_y_a[6*32 + 7]*stride_image + i_x_a[6*32 + 7]) < *(image_center + i_y_b[6*32 + 7]*stride_image + i_x_b[6*32 + 7])) << 7;
            f[6] |= (*(image_center + i_y_a[6*32 + 8]*stride_image + i_x_a[6*32 + 8]) < *(image_center + i_y_b[6*32 + 8]*stride_image + i_x_b[6*32 + 8])) << 8;
            f[6] |= (*(image_center + i_y_a[6*32 + 9]*stride_image + i_x_a[6*32 + 9]) < *(image_center + i_y_b[6*32 + 9]*stride_image + i_x_b[6*32 + 9])) << 9;
            f[6] |= (*(image_center + i_y_a[6*32 + 10]*stride_image + i_x_a[6*32 + 10]) < *(image_center + i_y_b[6*32 + 10]*stride_image + i_x_b[6*32 + 10])) << 10;
            f[6] |= (*(image_center + i_y_a[6*32 + 11]*stride_image + i_x_a[6*32 + 11]) < *(image_center + i_y_b[6*32 + 11]*stride_image + i_x_b[6*32 + 11])) << 11;
            f[6] |= (*(image_center + i_y_a[6*32 + 12]*stride_image + i_x_a[6*32 + 12]) < *(image_center + i_y_b[6*32 + 12]*stride_image + i_x_b[6*32 + 12])) << 12;
            f[6] |= (*(image_center + i_y_a[6*32 + 13]*stride_image + i_x_a[6*32 + 13]) < *(image_center + i_y_b[6*32 + 13]*stride_image + i_x_b[6*32 + 13])) << 13;
            f[6] |= (*(image_center + i_y_a[6*32 + 14]*stride_image + i_x_a[6*32 + 14]) < *(image_center + i_y_b[6*32 + 14]*stride_image + i_x_b[6*32 + 14])) << 14;
            f[6] |= (*(image_center + i_y_a[6*32 + 15]*stride_image + i_x_a[6*32 + 15]) < *(image_center + i_y_b[6*32 + 15]*stride_image + i_x_b[6*32 + 15])) << 15;
            f[6] |= (*(image_center + i_y_a[6*32 + 16]*stride_image + i_x_a[6*32 + 16]) < *(image_center + i_y_b[6*32 + 16]*stride_image + i_x_b[6*32 + 16])) << 16;
            f[6] |= (*(image_center + i_y_a[6*32 + 17]*stride_image + i_x_a[6*32 + 17]) < *(image_center + i_y_b[6*32 + 17]*stride_image + i_x_b[6*32 + 17])) << 17;
            f[6] |= (*(image_center + i_y_a[6*32 + 18]*stride_image + i_x_a[6*32 + 18]) < *(image_center + i_y_b[6*32 + 18]*stride_image + i_x_b[6*32 + 18])) << 18;
            f[6] |= (*(image_center + i_y_a[6*32 + 19]*stride_image + i_x_a[6*32 + 19]) < *(image_center + i_y_b[6*32 + 19]*stride_image + i_x_b[6*32 + 19])) << 19;
            f[6] |= (*(image_center + i_y_a[6*32 + 20]*stride_image + i_x_a[6*32 + 20]) < *(image_center + i_y_b[6*32 + 20]*stride_image + i_x_b[6*32 + 20])) << 20;
            f[6] |= (*(image_center + i_y_a[6*32 + 21]*stride_image + i_x_a[6*32 + 21]) < *(image_center + i_y_b[6*32 + 21]*stride_image + i_x_b[6*32 + 21])) << 21;
            f[6] |= (*(image_center + i_y_a[6*32 + 22]*stride_image + i_x_a[6*32 + 22]) < *(image_center + i_y_b[6*32 + 22]*stride_image + i_x_b[6*32 + 22])) << 22;
            f[6] |= (*(image_center + i_y_a[6*32 + 23]*stride_image + i_x_a[6*32 + 23]) < *(image_center + i_y_b[6*32 + 23]*stride_image + i_x_b[6*32 + 23])) << 23;
            f[6] |= (*(image_center + i_y_a[6*32 + 24]*stride_image + i_x_a[6*32 + 24]) < *(image_center + i_y_b[6*32 + 24]*stride_image + i_x_b[6*32 + 24])) << 24;
            f[6] |= (*(image_center + i_y_a[6*32 + 25]*stride_image + i_x_a[6*32 + 25]) < *(image_center + i_y_b[6*32 + 25]*stride_image + i_x_b[6*32 + 25])) << 25;
            f[6] |= (*(image_center + i_y_a[6*32 + 26]*stride_image + i_x_a[6*32 + 26]) < *(image_center + i_y_b[6*32 + 26]*stride_image + i_x_b[6*32 + 26])) << 26;
            f[6] |= (*(image_center + i_y_a[6*32 + 27]*stride_image + i_x_a[6*32 + 27]) < *(image_center + i_y_b[6*32 + 27]*stride_image + i_x_b[6*32 + 27])) << 27;
            f[6] |= (*(image_center + i_y_a[6*32 + 28]*stride_image + i_x_a[6*32 + 28]) < *(image_center + i_y_b[6*32 + 28]*stride_image + i_x_b[6*32 + 28])) << 28;
            f[6] |= (*(image_center + i_y_a[6*32 + 29]*stride_image + i_x_a[6*32 + 29]) < *(image_center + i_y_b[6*32 + 29]*stride_image + i_x_b[6*32 + 29])) << 29;
            f[6] |= (*(image_center + i_y_a[6*32 + 30]*stride_image + i_x_a[6*32 + 30]) < *(image_center + i_y_b[6*32 + 30]*stride_image + i_x_b[6*32 + 30])) << 30;
            f[6] |= (*(image_center + i_y_a[6*32 + 31]*stride_image + i_x_a[6*32 + 31]) < *(image_center + i_y_b[6*32 + 31]*stride_image + i_x_b[6*32 + 31])) << 31;
            
            f[7] = (*(image_center + i_y_a[7*32]*stride_image + i_x_a[7*32])
                        < *(image_center + i_y_b[7*32]*stride_image + i_x_b[7*32]));
            f[7] |= (*(image_center + i_y_a[7*32 + 1]*stride_image + i_x_a[7*32 + 1]) < *(image_center + i_y_b[7*32 + 1]*stride_image + i_x_b[7*32 + 1])) << 1;
            f[7] |= (*(image_center + i_y_a[7*32 + 2]*stride_image + i_x_a[7*32 + 2]) < *(image_center + i_y_b[7*32 + 2]*stride_image + i_x_b[7*32 + 2])) << 2;
            f[7] |= (*(image_center + i_y_a[7*32 + 3]*stride_image + i_x_a[7*32 + 3]) < *(image_center + i_y_b[7*32 + 3]*stride_image + i_x_b[7*32 + 3])) << 3;
            f[7] |= (*(image_center + i_y_a[7*32 + 4]*stride_image + i_x_a[7*32 + 4]) < *(image_center + i_y_b[7*32 + 4]*stride_image + i_x_b[7*32 + 4])) << 4;
            f[7] |= (*(image_center + i_y_a[7*32 + 5]*stride_image + i_x_a[7*32 + 5]) < *(image_center + i_y_b[7*32 + 5]*stride_image + i_x_b[7*32 + 5])) << 5;
            f[7] |= (*(image_center + i_y_a[7*32 + 6]*stride_image + i_x_a[7*32 + 6]) < *(image_center + i_y_b[7*32 + 6]*stride_image + i_x_b[7*32 + 6])) << 6;
            f[7] |= (*(image_center + i_y_a[7*32 + 7]*stride_image + i_x_a[7*32 + 7]) < *(image_center + i_y_b[7*32 + 7]*stride_image + i_x_b[7*32 + 7])) << 7;
            f[7] |= (*(image_center + i_y_a[7*32 + 8]*stride_image + i_x_a[7*32 + 8]) < *(image_center + i_y_b[7*32 + 8]*stride_image + i_x_b[7*32 + 8])) << 8;
            f[7] |= (*(image_center + i_y_a[7*32 + 9]*stride_image + i_x_a[7*32 + 9]) < *(image_center + i_y_b[7*32 + 9]*stride_image + i_x_b[7*32 + 9])) << 9;
            f[7] |= (*(image_center + i_y_a[7*32 + 10]*stride_image + i_x_a[7*32 + 10]) < *(image_center + i_y_b[7*32 + 10]*stride_image + i_x_b[7*32 + 10])) << 10;
            f[7] |= (*(image_center + i_y_a[7*32 + 11]*stride_image + i_x_a[7*32 + 11]) < *(image_center + i_y_b[7*32 + 11]*stride_image + i_x_b[7*32 + 11])) << 11;
            f[7] |= (*(image_center + i_y_a[7*32 + 12]*stride_image + i_x_a[7*32 + 12]) < *(image_center + i_y_b[7*32 + 12]*stride_image + i_x_b[7*32 + 12])) << 12;
            f[7] |= (*(image_center + i_y_a[7*32 + 13]*stride_image + i_x_a[7*32 + 13]) < *(image_center + i_y_b[7*32 + 13]*stride_image + i_x_b[7*32 + 13])) << 13;
            f[7] |= (*(image_center + i_y_a[7*32 + 14]*stride_image + i_x_a[7*32 + 14]) < *(image_center + i_y_b[7*32 + 14]*stride_image + i_x_b[7*32 + 14])) << 14;
            f[7] |= (*(image_center + i_y_a[7*32 + 15]*stride_image + i_x_a[7*32 + 15]) < *(image_center + i_y_b[7*32 + 15]*stride_image + i_x_b[7*32 + 15])) << 15;
            f[7] |= (*(image_center + i_y_a[7*32 + 16]*stride_image + i_x_a[7*32 + 16]) < *(image_center + i_y_b[7*32 + 16]*stride_image + i_x_b[7*32 + 16])) << 16;
            f[7] |= (*(image_center + i_y_a[7*32 + 17]*stride_image + i_x_a[7*32 + 17]) < *(image_center + i_y_b[7*32 + 17]*stride_image + i_x_b[7*32 + 17])) << 17;
            f[7] |= (*(image_center + i_y_a[7*32 + 18]*stride_image + i_x_a[7*32 + 18]) < *(image_center + i_y_b[7*32 + 18]*stride_image + i_x_b[7*32 + 18])) << 18;
            f[7] |= (*(image_center + i_y_a[7*32 + 19]*stride_image + i_x_a[7*32 + 19]) < *(image_center + i_y_b[7*32 + 19]*stride_image + i_x_b[7*32 + 19])) << 19;
            f[7] |= (*(image_center + i_y_a[7*32 + 20]*stride_image + i_x_a[7*32 + 20]) < *(image_center + i_y_b[7*32 + 20]*stride_image + i_x_b[7*32 + 20])) << 20;
            f[7] |= (*(image_center + i_y_a[7*32 + 21]*stride_image + i_x_a[7*32 + 21]) < *(image_center + i_y_b[7*32 + 21]*stride_image + i_x_b[7*32 + 21])) << 21;
            f[7] |= (*(image_center + i_y_a[7*32 + 22]*stride_image + i_x_a[7*32 + 22]) < *(image_center + i_y_b[7*32 + 22]*stride_image + i_x_b[7*32 + 22])) << 22;
            f[7] |= (*(image_center + i_y_a[7*32 + 23]*stride_image + i_x_a[7*32 + 23]) < *(image_center + i_y_b[7*32 + 23]*stride_image + i_x_b[7*32 + 23])) << 23;
            f[7] |= (*(image_center + i_y_a[7*32 + 24]*stride_image + i_x_a[7*32 + 24]) < *(image_center + i_y_b[7*32 + 24]*stride_image + i_x_b[7*32 + 24])) << 24;
            f[7] |= (*(image_center + i_y_a[7*32 + 25]*stride_image + i_x_a[7*32 + 25]) < *(image_center + i_y_b[7*32 + 25]*stride_image + i_x_b[7*32 + 25])) << 25;
            f[7] |= (*(image_center + i_y_a[7*32 + 26]*stride_image + i_x_a[7*32 + 26]) < *(image_center + i_y_b[7*32 + 26]*stride_image + i_x_b[7*32 + 26])) << 26;
            f[7] |= (*(image_center + i_y_a[7*32 + 27]*stride_image + i_x_a[7*32 + 27]) < *(image_center + i_y_b[7*32 + 27]*stride_image + i_x_b[7*32 + 27])) << 27;
            f[7] |= (*(image_center + i_y_a[7*32 + 28]*stride_image + i_x_a[7*32 + 28]) < *(image_center + i_y_b[7*32 + 28]*stride_image + i_x_b[7*32 + 28])) << 28;
            f[7] |= (*(image_center + i_y_a[7*32 + 29]*stride_image + i_x_a[7*32 + 29]) < *(image_center + i_y_b[7*32 + 29]*stride_image + i_x_b[7*32 + 29])) << 29;
            f[7] |= (*(image_center + i_y_a[7*32 + 30]*stride_image + i_x_a[7*32 + 30]) < *(image_center + i_y_b[7*32 + 30]*stride_image + i_x_b[7*32 + 30])) << 30;
            f[7] |= (*(image_center + i_y_a[7*32 + 31]*stride_image + i_x_a[7*32 + 31]) < *(image_center + i_y_b[7*32 + 31]*stride_image + i_x_b[7*32 + 31])) << 31;
            
            
            _mm_store_si128((__m128i*)(bd.memptr() + j*n_rows + 0*2),_mm_load_si128((const __m128i *)(f + 0*4)));
            _mm_store_si128((__m128i*)(bd.memptr() + j*n_rows + 1*2),_mm_load_si128((const __m128i *)(f + 1*4)));
            

        }
    }
}

    int n_rows;
    int n_cols;
    Matrix<int64_t> bd;
    static int diag_length_pattern; // <- maximal range of pattern box: 25/2 = 12, sqrt(12*12 + 12*12) = 17
    static int gaussian_bit_pattern_31_x_a[256] __attribute__((aligned(32)));
    static int gaussian_bit_pattern_31_y_a[256] __attribute__((aligned(32)));
    static int gaussian_bit_pattern_31_x_b[256] __attribute__((aligned(32)));
    static int gaussian_bit_pattern_31_y_b[256] __attribute__((aligned(32)));
};

#endif
