

#include "Types.h"
#include "BRIEF.h"
#include <cmath>
#include <iostream>
#include "xmmintrin.h"





#define GET_VALUE(i_pattern) \
               (xj = gaussian_bit_pattern_31[i_pattern]*cos_angle - gaussian_bit_pattern_31[i_pattern+1]*sin_angle, \
                yj = gaussian_bit_pattern_31[i_pattern]*sin_angle + gaussian_bit_pattern_31[i_pattern+1]*cos_angle, \
                ix = _mm_cvtss_si32(_mm_set_ss( xj )), \
                iy = _mm_cvtss_si32(_mm_set_ss( yj )), \
                *(image_src_center + iy*stride_image + ix*n_channels) )

void
BRIEF::rbrief(unsigned char *image_src, const int height_image, const int width_image, const int n_channels,
                const int stride_image, const int *x, const int *y, const float *angle, const int n_features, int64_t *bd,
                const int n_rows_bd) {
    for (int j = 0; j < n_features; j++) {
        if ((x[j] > diag_length_pattern) && x[j] < (width_image - diag_length_pattern)
            && (y[j] > diag_length_pattern) && y[j] < (height_image - diag_length_pattern)) {
            float cos_angle = std::cos(angle[j]);
            float sin_angle = std::sin(angle[j]);
            float xj, yj;
            int ix, iy;
            unsigned char *image_src_center = image_src + y[j] * stride_image + x[j] * n_channels;
            // N_DIM_BINARYDESCRIPTOR / SIZE_BITS_HAMING = 4
            for (int i = 0; i < N_DIM_BINARYDESCRIPTOR / SIZE_BITS_HAMING; i++) {
                int i_pattern = i * SIZE_BITS_HAMING*4;
                int64_t t0, t1, val;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val = t0 < t1;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 1;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 2;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 3;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 4;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 5;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 6;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 7;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 8;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 9;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 10;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 11;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 12;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 13;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 14;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 15;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 16;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 17;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 18;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 19;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 20;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 21;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 22;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 23;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 24;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 25;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 26;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 27;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 28;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 29;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 30;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 31;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 32;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 33;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 34;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 35;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 36;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 37;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 38;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 39;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 40;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 41;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 42;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 43;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 44;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 45;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 46;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 47;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 48;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 49;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 50;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 51;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 52;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 53;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 54;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 55;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 56;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 57;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 58;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 59;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 60;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 61;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 62;
                i_pattern += 4;
                t0 = GET_VALUE(i_pattern); t1 = GET_VALUE(i_pattern + 2);
                val |= ( int64_t ) (t0 < t1) << 63;
                ;
                bd[j*n_rows_bd + i] = val;
            }
        }
    }
}