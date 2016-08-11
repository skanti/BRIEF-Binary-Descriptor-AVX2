include(`./src/Unroll.m4')
#include "Types.h"
#include "BRIEF.h"
#include <cmath>
#include <iostream>
#include "xmmintrin.h"

define(intx_t,`ifelse(SIZE_BITS_HAMING,32, int32_t, int64_t)')
define(intx_suffix,`ifelse(SIZE_BITS_HAMING,32,,L)')
define(popcount_x,`ifelse(SIZE_BITS_HAMING,32,_popcnt32,_popcnt64)')

#define GET_VALUE(i_pattern) \
               (xj = gaussian_bit_pattern_31[i_pattern]*cos_angle - gaussian_bit_pattern_31[i_pattern+1]*sin_angle, \
                yj = gaussian_bit_pattern_31[i_pattern]*sin_angle + gaussian_bit_pattern_31[i_pattern+1]*cos_angle, \
                ix = _mm_cvtss_si32(_mm_set_ss( xj )), \
                iy = _mm_cvtss_si32(_mm_set_ss( yj )), \
                *(image_src_center + iy*stride_image + ix*n_channels))

void
BRIEF::rbrief(unsigned char *image_src, const int height_image, const int width_image, const int n_channels,
                const int stride_image, const int *x, const int *y, const float *angle, const int n_features, intx_t *bd,
                const int n_rows_bd) {
    for (int j = 0; j < n_features; j++) {
        if ((x[j] > diag_length_pattern) && x[j] < (width_image - diag_length_pattern)
            && (y[j] > diag_length_pattern) && y[j] < (height_image - diag_length_pattern)) {
            float cos_angle = std::cos(angle[j]);
            float sin_angle = std::sin(angle[j]);
            float xj, yj;
            int ix, iy;
            unsigned char *image_src_center = image_src + y[j] * stride_image + x[j] * n_channels;
            // `N_DIM_BINARYDESCRIPTOR' / `SIZE_BITS_HAMING' = eval( N_DIM_BINARYDESCRIPTOR / SIZE_BITS_HAMING)
            for (int i = 0; i < `N_DIM_BINARYDESCRIPTOR' / `SIZE_BITS_HAMING'; i++) {
                int i_pattern = i * `SIZE_BITS_HAMING';
                intx_t t0, t1, val;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val = t0 < t1;
                forloop(k,1,eval(SIZE_BITS_HAMING-1),
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( intx_t ) (t0 < t1) << `k';
                );
                bd[j*n_rows_bd + i] = val;
            }
        }
    }
}