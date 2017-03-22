include(`./src/Unroll.m4')#pragma once

#include "Types.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include "immintrin.h"

#ifndef N_DIM_BINARYDESCRIPTOR
#define N_DIM_BINARYDESCRIPTOR 256
#else
assert(`N_DIM_BINARYDESCRIPTOR' == 256)
#endif

#ifndef SIZE_BITS_HAMING
#define SIZE_BITS_HAMING 64
#else
assert(`SIZE_BITS_HAMING' == 64)
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
            // `N_DIM_BINARYDESCRIPTOR' / `SIZE_BITS_HAMING' = eval( N_DIM_BINARYDESCRIPTOR / SIZE_BITS_HAMING)
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
            forloop(k,0,7,
            f[k] = (*(image_center + i_y_a[k*32]*stride_image + i_x_a[k*32])
                        < *(image_center + i_y_b[k*32]*stride_image + i_x_b[k*32]));
            `forloop(l,1,31,
            f[k] |= (*(image_center + i_y_a[k*32 + l]*stride_image + i_x_a[k*32 + l]) < *(image_center + i_y_b[k*32 + l]*stride_image + i_x_b[k*32 + l])) << l;
            )'
            )
            forloop(k,0,1,
            _mm_store_si128((__m128i*)(bd.memptr() + j*n_rows + k*2),_mm_load_si128((const __m128i *)(f + k*4)));
            )

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
