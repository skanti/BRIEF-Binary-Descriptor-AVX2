#ifndef BRIEF_H
#define BRIEF_H

#include "Types.h"
#include <cassert>

#ifndef N_DIM_BINARYDESCRIPTOR
#define N_DIM_BINARYDESCRIPTOR 256
#else
assert(N_DIM_BINARYDESCRIPTOR == 256)
#endif

#ifndef SIZE_BITS_HAMING
#define SIZE_BITS_HAMING 64
#else
assert(SIZE_BITS_HAMING == 64)
#endif

class BRIEF {
public:
    static void rbrief(unsigned char *image_src, const int height_image, const int width_image, const int n_channels,
                       const int stride_image, const int *x, const int *y, const float *angle, const int n_features,
                       int64_t *bd, const int n_rows_bd);

    static int diag_length_pattern; // <- maximal range of pattern box: 25/2 = 12, sqrt(12*12 + 12*12) = 17
    static int gaussian_bit_pattern_31_x_a[256] __attribute__((aligned(32)));
    static int gaussian_bit_pattern_31_y_a[256] __attribute__((aligned(32)));
    static int gaussian_bit_pattern_31_x_b[256] __attribute__((aligned(32)));
    static int gaussian_bit_pattern_31_y_b[256] __attribute__((aligned(32)));
};

#endif
