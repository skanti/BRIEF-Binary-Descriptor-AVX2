

#include "Types.h"
#include "BRIEF.h"
#include <cmath>
#include <iostream>
#include "xmmintrin.h"
#include "emmintrin.h"





#define GET_VALUE(k, x_f,y_f, x_i, y_i) \
                ((x_f) = gaussian_bit_pattern_31[k]*cos_angle - gaussian_bit_pattern_31[k+1]*sin_angle, \
                (y_f) = gaussian_bit_pattern_31[k]*sin_angle + gaussian_bit_pattern_31[k+1]*cos_angle, \
                (x_i) = _mm_cvtss_si32(_mm_set_ss( (x_f) )), \
                (y_i) = _mm_cvtss_si32(_mm_set_ss( (y_f) )), \
                *(image_src_center + y_i*stride_image + x_i*n_channels))

void
BRIEF::rbrief(unsigned char *image_src, const int height_image, const int width_image, const int n_channels,
                const int stride_image, const int *x, const int *y, const float *angle, const int n_features, int64_t *bd,
                const int n_rows_bd) {
    for (int j = 0; j < n_features; j++) {
        if ((x[j] > diag_length_pattern) && x[j] < (width_image - diag_length_pattern)
            && (y[j] > diag_length_pattern) && y[j] < (height_image - diag_length_pattern)) {
            float cos_angle = std::cos(angle[j]);
            float sin_angle = std::sin(angle[j]);
            unsigned char *image_src_center = image_src + y[j] * stride_image + x[j] * n_channels;
            // N_DIM_BINARYDESCRIPTOR / SIZE_BITS_HAMING = 4
            unsigned char a[256] __attribute__((aligned(32)));
            unsigned char b[256] __attribute__((aligned(32)));
            float x_af, x_bf, y_af, y_bf;
            int x_a, x_b, y_a, y_b;
            a[0] = GET_VALUE(4*0,x_af,y_af,x_a,y_a); b[0] = GET_VALUE(0*4 + 2,x_bf,y_bf,x_b,y_b);
            a[1] = GET_VALUE(4*1,x_af,y_af,x_a,y_a); b[1] = GET_VALUE(1*4 + 2,x_bf,y_bf,x_b,y_b);
            a[2] = GET_VALUE(4*2,x_af,y_af,x_a,y_a); b[2] = GET_VALUE(2*4 + 2,x_bf,y_bf,x_b,y_b);
            a[3] = GET_VALUE(4*3,x_af,y_af,x_a,y_a); b[3] = GET_VALUE(3*4 + 2,x_bf,y_bf,x_b,y_b);
            a[4] = GET_VALUE(4*4,x_af,y_af,x_a,y_a); b[4] = GET_VALUE(4*4 + 2,x_bf,y_bf,x_b,y_b);
            a[5] = GET_VALUE(4*5,x_af,y_af,x_a,y_a); b[5] = GET_VALUE(5*4 + 2,x_bf,y_bf,x_b,y_b);
            a[6] = GET_VALUE(4*6,x_af,y_af,x_a,y_a); b[6] = GET_VALUE(6*4 + 2,x_bf,y_bf,x_b,y_b);
            a[7] = GET_VALUE(4*7,x_af,y_af,x_a,y_a); b[7] = GET_VALUE(7*4 + 2,x_bf,y_bf,x_b,y_b);
            a[8] = GET_VALUE(4*8,x_af,y_af,x_a,y_a); b[8] = GET_VALUE(8*4 + 2,x_bf,y_bf,x_b,y_b);
            a[9] = GET_VALUE(4*9,x_af,y_af,x_a,y_a); b[9] = GET_VALUE(9*4 + 2,x_bf,y_bf,x_b,y_b);
            a[10] = GET_VALUE(4*10,x_af,y_af,x_a,y_a); b[10] = GET_VALUE(10*4 + 2,x_bf,y_bf,x_b,y_b);
            a[11] = GET_VALUE(4*11,x_af,y_af,x_a,y_a); b[11] = GET_VALUE(11*4 + 2,x_bf,y_bf,x_b,y_b);
            a[12] = GET_VALUE(4*12,x_af,y_af,x_a,y_a); b[12] = GET_VALUE(12*4 + 2,x_bf,y_bf,x_b,y_b);
            a[13] = GET_VALUE(4*13,x_af,y_af,x_a,y_a); b[13] = GET_VALUE(13*4 + 2,x_bf,y_bf,x_b,y_b);
            a[14] = GET_VALUE(4*14,x_af,y_af,x_a,y_a); b[14] = GET_VALUE(14*4 + 2,x_bf,y_bf,x_b,y_b);
            a[15] = GET_VALUE(4*15,x_af,y_af,x_a,y_a); b[15] = GET_VALUE(15*4 + 2,x_bf,y_bf,x_b,y_b);
            a[16] = GET_VALUE(4*16,x_af,y_af,x_a,y_a); b[16] = GET_VALUE(16*4 + 2,x_bf,y_bf,x_b,y_b);
            a[17] = GET_VALUE(4*17,x_af,y_af,x_a,y_a); b[17] = GET_VALUE(17*4 + 2,x_bf,y_bf,x_b,y_b);
            a[18] = GET_VALUE(4*18,x_af,y_af,x_a,y_a); b[18] = GET_VALUE(18*4 + 2,x_bf,y_bf,x_b,y_b);
            a[19] = GET_VALUE(4*19,x_af,y_af,x_a,y_a); b[19] = GET_VALUE(19*4 + 2,x_bf,y_bf,x_b,y_b);
            a[20] = GET_VALUE(4*20,x_af,y_af,x_a,y_a); b[20] = GET_VALUE(20*4 + 2,x_bf,y_bf,x_b,y_b);
            a[21] = GET_VALUE(4*21,x_af,y_af,x_a,y_a); b[21] = GET_VALUE(21*4 + 2,x_bf,y_bf,x_b,y_b);
            a[22] = GET_VALUE(4*22,x_af,y_af,x_a,y_a); b[22] = GET_VALUE(22*4 + 2,x_bf,y_bf,x_b,y_b);
            a[23] = GET_VALUE(4*23,x_af,y_af,x_a,y_a); b[23] = GET_VALUE(23*4 + 2,x_bf,y_bf,x_b,y_b);
            a[24] = GET_VALUE(4*24,x_af,y_af,x_a,y_a); b[24] = GET_VALUE(24*4 + 2,x_bf,y_bf,x_b,y_b);
            a[25] = GET_VALUE(4*25,x_af,y_af,x_a,y_a); b[25] = GET_VALUE(25*4 + 2,x_bf,y_bf,x_b,y_b);
            a[26] = GET_VALUE(4*26,x_af,y_af,x_a,y_a); b[26] = GET_VALUE(26*4 + 2,x_bf,y_bf,x_b,y_b);
            a[27] = GET_VALUE(4*27,x_af,y_af,x_a,y_a); b[27] = GET_VALUE(27*4 + 2,x_bf,y_bf,x_b,y_b);
            a[28] = GET_VALUE(4*28,x_af,y_af,x_a,y_a); b[28] = GET_VALUE(28*4 + 2,x_bf,y_bf,x_b,y_b);
            a[29] = GET_VALUE(4*29,x_af,y_af,x_a,y_a); b[29] = GET_VALUE(29*4 + 2,x_bf,y_bf,x_b,y_b);
            a[30] = GET_VALUE(4*30,x_af,y_af,x_a,y_a); b[30] = GET_VALUE(30*4 + 2,x_bf,y_bf,x_b,y_b);
            a[31] = GET_VALUE(4*31,x_af,y_af,x_a,y_a); b[31] = GET_VALUE(31*4 + 2,x_bf,y_bf,x_b,y_b);
            a[32] = GET_VALUE(4*32,x_af,y_af,x_a,y_a); b[32] = GET_VALUE(32*4 + 2,x_bf,y_bf,x_b,y_b);
            a[33] = GET_VALUE(4*33,x_af,y_af,x_a,y_a); b[33] = GET_VALUE(33*4 + 2,x_bf,y_bf,x_b,y_b);
            a[34] = GET_VALUE(4*34,x_af,y_af,x_a,y_a); b[34] = GET_VALUE(34*4 + 2,x_bf,y_bf,x_b,y_b);
            a[35] = GET_VALUE(4*35,x_af,y_af,x_a,y_a); b[35] = GET_VALUE(35*4 + 2,x_bf,y_bf,x_b,y_b);
            a[36] = GET_VALUE(4*36,x_af,y_af,x_a,y_a); b[36] = GET_VALUE(36*4 + 2,x_bf,y_bf,x_b,y_b);
            a[37] = GET_VALUE(4*37,x_af,y_af,x_a,y_a); b[37] = GET_VALUE(37*4 + 2,x_bf,y_bf,x_b,y_b);
            a[38] = GET_VALUE(4*38,x_af,y_af,x_a,y_a); b[38] = GET_VALUE(38*4 + 2,x_bf,y_bf,x_b,y_b);
            a[39] = GET_VALUE(4*39,x_af,y_af,x_a,y_a); b[39] = GET_VALUE(39*4 + 2,x_bf,y_bf,x_b,y_b);
            a[40] = GET_VALUE(4*40,x_af,y_af,x_a,y_a); b[40] = GET_VALUE(40*4 + 2,x_bf,y_bf,x_b,y_b);
            a[41] = GET_VALUE(4*41,x_af,y_af,x_a,y_a); b[41] = GET_VALUE(41*4 + 2,x_bf,y_bf,x_b,y_b);
            a[42] = GET_VALUE(4*42,x_af,y_af,x_a,y_a); b[42] = GET_VALUE(42*4 + 2,x_bf,y_bf,x_b,y_b);
            a[43] = GET_VALUE(4*43,x_af,y_af,x_a,y_a); b[43] = GET_VALUE(43*4 + 2,x_bf,y_bf,x_b,y_b);
            a[44] = GET_VALUE(4*44,x_af,y_af,x_a,y_a); b[44] = GET_VALUE(44*4 + 2,x_bf,y_bf,x_b,y_b);
            a[45] = GET_VALUE(4*45,x_af,y_af,x_a,y_a); b[45] = GET_VALUE(45*4 + 2,x_bf,y_bf,x_b,y_b);
            a[46] = GET_VALUE(4*46,x_af,y_af,x_a,y_a); b[46] = GET_VALUE(46*4 + 2,x_bf,y_bf,x_b,y_b);
            a[47] = GET_VALUE(4*47,x_af,y_af,x_a,y_a); b[47] = GET_VALUE(47*4 + 2,x_bf,y_bf,x_b,y_b);
            a[48] = GET_VALUE(4*48,x_af,y_af,x_a,y_a); b[48] = GET_VALUE(48*4 + 2,x_bf,y_bf,x_b,y_b);
            a[49] = GET_VALUE(4*49,x_af,y_af,x_a,y_a); b[49] = GET_VALUE(49*4 + 2,x_bf,y_bf,x_b,y_b);
            a[50] = GET_VALUE(4*50,x_af,y_af,x_a,y_a); b[50] = GET_VALUE(50*4 + 2,x_bf,y_bf,x_b,y_b);
            a[51] = GET_VALUE(4*51,x_af,y_af,x_a,y_a); b[51] = GET_VALUE(51*4 + 2,x_bf,y_bf,x_b,y_b);
            a[52] = GET_VALUE(4*52,x_af,y_af,x_a,y_a); b[52] = GET_VALUE(52*4 + 2,x_bf,y_bf,x_b,y_b);
            a[53] = GET_VALUE(4*53,x_af,y_af,x_a,y_a); b[53] = GET_VALUE(53*4 + 2,x_bf,y_bf,x_b,y_b);
            a[54] = GET_VALUE(4*54,x_af,y_af,x_a,y_a); b[54] = GET_VALUE(54*4 + 2,x_bf,y_bf,x_b,y_b);
            a[55] = GET_VALUE(4*55,x_af,y_af,x_a,y_a); b[55] = GET_VALUE(55*4 + 2,x_bf,y_bf,x_b,y_b);
            a[56] = GET_VALUE(4*56,x_af,y_af,x_a,y_a); b[56] = GET_VALUE(56*4 + 2,x_bf,y_bf,x_b,y_b);
            a[57] = GET_VALUE(4*57,x_af,y_af,x_a,y_a); b[57] = GET_VALUE(57*4 + 2,x_bf,y_bf,x_b,y_b);
            a[58] = GET_VALUE(4*58,x_af,y_af,x_a,y_a); b[58] = GET_VALUE(58*4 + 2,x_bf,y_bf,x_b,y_b);
            a[59] = GET_VALUE(4*59,x_af,y_af,x_a,y_a); b[59] = GET_VALUE(59*4 + 2,x_bf,y_bf,x_b,y_b);
            a[60] = GET_VALUE(4*60,x_af,y_af,x_a,y_a); b[60] = GET_VALUE(60*4 + 2,x_bf,y_bf,x_b,y_b);
            a[61] = GET_VALUE(4*61,x_af,y_af,x_a,y_a); b[61] = GET_VALUE(61*4 + 2,x_bf,y_bf,x_b,y_b);
            a[62] = GET_VALUE(4*62,x_af,y_af,x_a,y_a); b[62] = GET_VALUE(62*4 + 2,x_bf,y_bf,x_b,y_b);
            a[63] = GET_VALUE(4*63,x_af,y_af,x_a,y_a); b[63] = GET_VALUE(63*4 + 2,x_bf,y_bf,x_b,y_b);
            a[64] = GET_VALUE(4*64,x_af,y_af,x_a,y_a); b[64] = GET_VALUE(64*4 + 2,x_bf,y_bf,x_b,y_b);
            a[65] = GET_VALUE(4*65,x_af,y_af,x_a,y_a); b[65] = GET_VALUE(65*4 + 2,x_bf,y_bf,x_b,y_b);
            a[66] = GET_VALUE(4*66,x_af,y_af,x_a,y_a); b[66] = GET_VALUE(66*4 + 2,x_bf,y_bf,x_b,y_b);
            a[67] = GET_VALUE(4*67,x_af,y_af,x_a,y_a); b[67] = GET_VALUE(67*4 + 2,x_bf,y_bf,x_b,y_b);
            a[68] = GET_VALUE(4*68,x_af,y_af,x_a,y_a); b[68] = GET_VALUE(68*4 + 2,x_bf,y_bf,x_b,y_b);
            a[69] = GET_VALUE(4*69,x_af,y_af,x_a,y_a); b[69] = GET_VALUE(69*4 + 2,x_bf,y_bf,x_b,y_b);
            a[70] = GET_VALUE(4*70,x_af,y_af,x_a,y_a); b[70] = GET_VALUE(70*4 + 2,x_bf,y_bf,x_b,y_b);
            a[71] = GET_VALUE(4*71,x_af,y_af,x_a,y_a); b[71] = GET_VALUE(71*4 + 2,x_bf,y_bf,x_b,y_b);
            a[72] = GET_VALUE(4*72,x_af,y_af,x_a,y_a); b[72] = GET_VALUE(72*4 + 2,x_bf,y_bf,x_b,y_b);
            a[73] = GET_VALUE(4*73,x_af,y_af,x_a,y_a); b[73] = GET_VALUE(73*4 + 2,x_bf,y_bf,x_b,y_b);
            a[74] = GET_VALUE(4*74,x_af,y_af,x_a,y_a); b[74] = GET_VALUE(74*4 + 2,x_bf,y_bf,x_b,y_b);
            a[75] = GET_VALUE(4*75,x_af,y_af,x_a,y_a); b[75] = GET_VALUE(75*4 + 2,x_bf,y_bf,x_b,y_b);
            a[76] = GET_VALUE(4*76,x_af,y_af,x_a,y_a); b[76] = GET_VALUE(76*4 + 2,x_bf,y_bf,x_b,y_b);
            a[77] = GET_VALUE(4*77,x_af,y_af,x_a,y_a); b[77] = GET_VALUE(77*4 + 2,x_bf,y_bf,x_b,y_b);
            a[78] = GET_VALUE(4*78,x_af,y_af,x_a,y_a); b[78] = GET_VALUE(78*4 + 2,x_bf,y_bf,x_b,y_b);
            a[79] = GET_VALUE(4*79,x_af,y_af,x_a,y_a); b[79] = GET_VALUE(79*4 + 2,x_bf,y_bf,x_b,y_b);
            a[80] = GET_VALUE(4*80,x_af,y_af,x_a,y_a); b[80] = GET_VALUE(80*4 + 2,x_bf,y_bf,x_b,y_b);
            a[81] = GET_VALUE(4*81,x_af,y_af,x_a,y_a); b[81] = GET_VALUE(81*4 + 2,x_bf,y_bf,x_b,y_b);
            a[82] = GET_VALUE(4*82,x_af,y_af,x_a,y_a); b[82] = GET_VALUE(82*4 + 2,x_bf,y_bf,x_b,y_b);
            a[83] = GET_VALUE(4*83,x_af,y_af,x_a,y_a); b[83] = GET_VALUE(83*4 + 2,x_bf,y_bf,x_b,y_b);
            a[84] = GET_VALUE(4*84,x_af,y_af,x_a,y_a); b[84] = GET_VALUE(84*4 + 2,x_bf,y_bf,x_b,y_b);
            a[85] = GET_VALUE(4*85,x_af,y_af,x_a,y_a); b[85] = GET_VALUE(85*4 + 2,x_bf,y_bf,x_b,y_b);
            a[86] = GET_VALUE(4*86,x_af,y_af,x_a,y_a); b[86] = GET_VALUE(86*4 + 2,x_bf,y_bf,x_b,y_b);
            a[87] = GET_VALUE(4*87,x_af,y_af,x_a,y_a); b[87] = GET_VALUE(87*4 + 2,x_bf,y_bf,x_b,y_b);
            a[88] = GET_VALUE(4*88,x_af,y_af,x_a,y_a); b[88] = GET_VALUE(88*4 + 2,x_bf,y_bf,x_b,y_b);
            a[89] = GET_VALUE(4*89,x_af,y_af,x_a,y_a); b[89] = GET_VALUE(89*4 + 2,x_bf,y_bf,x_b,y_b);
            a[90] = GET_VALUE(4*90,x_af,y_af,x_a,y_a); b[90] = GET_VALUE(90*4 + 2,x_bf,y_bf,x_b,y_b);
            a[91] = GET_VALUE(4*91,x_af,y_af,x_a,y_a); b[91] = GET_VALUE(91*4 + 2,x_bf,y_bf,x_b,y_b);
            a[92] = GET_VALUE(4*92,x_af,y_af,x_a,y_a); b[92] = GET_VALUE(92*4 + 2,x_bf,y_bf,x_b,y_b);
            a[93] = GET_VALUE(4*93,x_af,y_af,x_a,y_a); b[93] = GET_VALUE(93*4 + 2,x_bf,y_bf,x_b,y_b);
            a[94] = GET_VALUE(4*94,x_af,y_af,x_a,y_a); b[94] = GET_VALUE(94*4 + 2,x_bf,y_bf,x_b,y_b);
            a[95] = GET_VALUE(4*95,x_af,y_af,x_a,y_a); b[95] = GET_VALUE(95*4 + 2,x_bf,y_bf,x_b,y_b);
            a[96] = GET_VALUE(4*96,x_af,y_af,x_a,y_a); b[96] = GET_VALUE(96*4 + 2,x_bf,y_bf,x_b,y_b);
            a[97] = GET_VALUE(4*97,x_af,y_af,x_a,y_a); b[97] = GET_VALUE(97*4 + 2,x_bf,y_bf,x_b,y_b);
            a[98] = GET_VALUE(4*98,x_af,y_af,x_a,y_a); b[98] = GET_VALUE(98*4 + 2,x_bf,y_bf,x_b,y_b);
            a[99] = GET_VALUE(4*99,x_af,y_af,x_a,y_a); b[99] = GET_VALUE(99*4 + 2,x_bf,y_bf,x_b,y_b);
            a[100] = GET_VALUE(4*100,x_af,y_af,x_a,y_a); b[100] = GET_VALUE(100*4 + 2,x_bf,y_bf,x_b,y_b);
            a[101] = GET_VALUE(4*101,x_af,y_af,x_a,y_a); b[101] = GET_VALUE(101*4 + 2,x_bf,y_bf,x_b,y_b);
            a[102] = GET_VALUE(4*102,x_af,y_af,x_a,y_a); b[102] = GET_VALUE(102*4 + 2,x_bf,y_bf,x_b,y_b);
            a[103] = GET_VALUE(4*103,x_af,y_af,x_a,y_a); b[103] = GET_VALUE(103*4 + 2,x_bf,y_bf,x_b,y_b);
            a[104] = GET_VALUE(4*104,x_af,y_af,x_a,y_a); b[104] = GET_VALUE(104*4 + 2,x_bf,y_bf,x_b,y_b);
            a[105] = GET_VALUE(4*105,x_af,y_af,x_a,y_a); b[105] = GET_VALUE(105*4 + 2,x_bf,y_bf,x_b,y_b);
            a[106] = GET_VALUE(4*106,x_af,y_af,x_a,y_a); b[106] = GET_VALUE(106*4 + 2,x_bf,y_bf,x_b,y_b);
            a[107] = GET_VALUE(4*107,x_af,y_af,x_a,y_a); b[107] = GET_VALUE(107*4 + 2,x_bf,y_bf,x_b,y_b);
            a[108] = GET_VALUE(4*108,x_af,y_af,x_a,y_a); b[108] = GET_VALUE(108*4 + 2,x_bf,y_bf,x_b,y_b);
            a[109] = GET_VALUE(4*109,x_af,y_af,x_a,y_a); b[109] = GET_VALUE(109*4 + 2,x_bf,y_bf,x_b,y_b);
            a[110] = GET_VALUE(4*110,x_af,y_af,x_a,y_a); b[110] = GET_VALUE(110*4 + 2,x_bf,y_bf,x_b,y_b);
            a[111] = GET_VALUE(4*111,x_af,y_af,x_a,y_a); b[111] = GET_VALUE(111*4 + 2,x_bf,y_bf,x_b,y_b);
            a[112] = GET_VALUE(4*112,x_af,y_af,x_a,y_a); b[112] = GET_VALUE(112*4 + 2,x_bf,y_bf,x_b,y_b);
            a[113] = GET_VALUE(4*113,x_af,y_af,x_a,y_a); b[113] = GET_VALUE(113*4 + 2,x_bf,y_bf,x_b,y_b);
            a[114] = GET_VALUE(4*114,x_af,y_af,x_a,y_a); b[114] = GET_VALUE(114*4 + 2,x_bf,y_bf,x_b,y_b);
            a[115] = GET_VALUE(4*115,x_af,y_af,x_a,y_a); b[115] = GET_VALUE(115*4 + 2,x_bf,y_bf,x_b,y_b);
            a[116] = GET_VALUE(4*116,x_af,y_af,x_a,y_a); b[116] = GET_VALUE(116*4 + 2,x_bf,y_bf,x_b,y_b);
            a[117] = GET_VALUE(4*117,x_af,y_af,x_a,y_a); b[117] = GET_VALUE(117*4 + 2,x_bf,y_bf,x_b,y_b);
            a[118] = GET_VALUE(4*118,x_af,y_af,x_a,y_a); b[118] = GET_VALUE(118*4 + 2,x_bf,y_bf,x_b,y_b);
            a[119] = GET_VALUE(4*119,x_af,y_af,x_a,y_a); b[119] = GET_VALUE(119*4 + 2,x_bf,y_bf,x_b,y_b);
            a[120] = GET_VALUE(4*120,x_af,y_af,x_a,y_a); b[120] = GET_VALUE(120*4 + 2,x_bf,y_bf,x_b,y_b);
            a[121] = GET_VALUE(4*121,x_af,y_af,x_a,y_a); b[121] = GET_VALUE(121*4 + 2,x_bf,y_bf,x_b,y_b);
            a[122] = GET_VALUE(4*122,x_af,y_af,x_a,y_a); b[122] = GET_VALUE(122*4 + 2,x_bf,y_bf,x_b,y_b);
            a[123] = GET_VALUE(4*123,x_af,y_af,x_a,y_a); b[123] = GET_VALUE(123*4 + 2,x_bf,y_bf,x_b,y_b);
            a[124] = GET_VALUE(4*124,x_af,y_af,x_a,y_a); b[124] = GET_VALUE(124*4 + 2,x_bf,y_bf,x_b,y_b);
            a[125] = GET_VALUE(4*125,x_af,y_af,x_a,y_a); b[125] = GET_VALUE(125*4 + 2,x_bf,y_bf,x_b,y_b);
            a[126] = GET_VALUE(4*126,x_af,y_af,x_a,y_a); b[126] = GET_VALUE(126*4 + 2,x_bf,y_bf,x_b,y_b);
            a[127] = GET_VALUE(4*127,x_af,y_af,x_a,y_a); b[127] = GET_VALUE(127*4 + 2,x_bf,y_bf,x_b,y_b);
            a[128] = GET_VALUE(4*128,x_af,y_af,x_a,y_a); b[128] = GET_VALUE(128*4 + 2,x_bf,y_bf,x_b,y_b);
            a[129] = GET_VALUE(4*129,x_af,y_af,x_a,y_a); b[129] = GET_VALUE(129*4 + 2,x_bf,y_bf,x_b,y_b);
            a[130] = GET_VALUE(4*130,x_af,y_af,x_a,y_a); b[130] = GET_VALUE(130*4 + 2,x_bf,y_bf,x_b,y_b);
            a[131] = GET_VALUE(4*131,x_af,y_af,x_a,y_a); b[131] = GET_VALUE(131*4 + 2,x_bf,y_bf,x_b,y_b);
            a[132] = GET_VALUE(4*132,x_af,y_af,x_a,y_a); b[132] = GET_VALUE(132*4 + 2,x_bf,y_bf,x_b,y_b);
            a[133] = GET_VALUE(4*133,x_af,y_af,x_a,y_a); b[133] = GET_VALUE(133*4 + 2,x_bf,y_bf,x_b,y_b);
            a[134] = GET_VALUE(4*134,x_af,y_af,x_a,y_a); b[134] = GET_VALUE(134*4 + 2,x_bf,y_bf,x_b,y_b);
            a[135] = GET_VALUE(4*135,x_af,y_af,x_a,y_a); b[135] = GET_VALUE(135*4 + 2,x_bf,y_bf,x_b,y_b);
            a[136] = GET_VALUE(4*136,x_af,y_af,x_a,y_a); b[136] = GET_VALUE(136*4 + 2,x_bf,y_bf,x_b,y_b);
            a[137] = GET_VALUE(4*137,x_af,y_af,x_a,y_a); b[137] = GET_VALUE(137*4 + 2,x_bf,y_bf,x_b,y_b);
            a[138] = GET_VALUE(4*138,x_af,y_af,x_a,y_a); b[138] = GET_VALUE(138*4 + 2,x_bf,y_bf,x_b,y_b);
            a[139] = GET_VALUE(4*139,x_af,y_af,x_a,y_a); b[139] = GET_VALUE(139*4 + 2,x_bf,y_bf,x_b,y_b);
            a[140] = GET_VALUE(4*140,x_af,y_af,x_a,y_a); b[140] = GET_VALUE(140*4 + 2,x_bf,y_bf,x_b,y_b);
            a[141] = GET_VALUE(4*141,x_af,y_af,x_a,y_a); b[141] = GET_VALUE(141*4 + 2,x_bf,y_bf,x_b,y_b);
            a[142] = GET_VALUE(4*142,x_af,y_af,x_a,y_a); b[142] = GET_VALUE(142*4 + 2,x_bf,y_bf,x_b,y_b);
            a[143] = GET_VALUE(4*143,x_af,y_af,x_a,y_a); b[143] = GET_VALUE(143*4 + 2,x_bf,y_bf,x_b,y_b);
            a[144] = GET_VALUE(4*144,x_af,y_af,x_a,y_a); b[144] = GET_VALUE(144*4 + 2,x_bf,y_bf,x_b,y_b);
            a[145] = GET_VALUE(4*145,x_af,y_af,x_a,y_a); b[145] = GET_VALUE(145*4 + 2,x_bf,y_bf,x_b,y_b);
            a[146] = GET_VALUE(4*146,x_af,y_af,x_a,y_a); b[146] = GET_VALUE(146*4 + 2,x_bf,y_bf,x_b,y_b);
            a[147] = GET_VALUE(4*147,x_af,y_af,x_a,y_a); b[147] = GET_VALUE(147*4 + 2,x_bf,y_bf,x_b,y_b);
            a[148] = GET_VALUE(4*148,x_af,y_af,x_a,y_a); b[148] = GET_VALUE(148*4 + 2,x_bf,y_bf,x_b,y_b);
            a[149] = GET_VALUE(4*149,x_af,y_af,x_a,y_a); b[149] = GET_VALUE(149*4 + 2,x_bf,y_bf,x_b,y_b);
            a[150] = GET_VALUE(4*150,x_af,y_af,x_a,y_a); b[150] = GET_VALUE(150*4 + 2,x_bf,y_bf,x_b,y_b);
            a[151] = GET_VALUE(4*151,x_af,y_af,x_a,y_a); b[151] = GET_VALUE(151*4 + 2,x_bf,y_bf,x_b,y_b);
            a[152] = GET_VALUE(4*152,x_af,y_af,x_a,y_a); b[152] = GET_VALUE(152*4 + 2,x_bf,y_bf,x_b,y_b);
            a[153] = GET_VALUE(4*153,x_af,y_af,x_a,y_a); b[153] = GET_VALUE(153*4 + 2,x_bf,y_bf,x_b,y_b);
            a[154] = GET_VALUE(4*154,x_af,y_af,x_a,y_a); b[154] = GET_VALUE(154*4 + 2,x_bf,y_bf,x_b,y_b);
            a[155] = GET_VALUE(4*155,x_af,y_af,x_a,y_a); b[155] = GET_VALUE(155*4 + 2,x_bf,y_bf,x_b,y_b);
            a[156] = GET_VALUE(4*156,x_af,y_af,x_a,y_a); b[156] = GET_VALUE(156*4 + 2,x_bf,y_bf,x_b,y_b);
            a[157] = GET_VALUE(4*157,x_af,y_af,x_a,y_a); b[157] = GET_VALUE(157*4 + 2,x_bf,y_bf,x_b,y_b);
            a[158] = GET_VALUE(4*158,x_af,y_af,x_a,y_a); b[158] = GET_VALUE(158*4 + 2,x_bf,y_bf,x_b,y_b);
            a[159] = GET_VALUE(4*159,x_af,y_af,x_a,y_a); b[159] = GET_VALUE(159*4 + 2,x_bf,y_bf,x_b,y_b);
            a[160] = GET_VALUE(4*160,x_af,y_af,x_a,y_a); b[160] = GET_VALUE(160*4 + 2,x_bf,y_bf,x_b,y_b);
            a[161] = GET_VALUE(4*161,x_af,y_af,x_a,y_a); b[161] = GET_VALUE(161*4 + 2,x_bf,y_bf,x_b,y_b);
            a[162] = GET_VALUE(4*162,x_af,y_af,x_a,y_a); b[162] = GET_VALUE(162*4 + 2,x_bf,y_bf,x_b,y_b);
            a[163] = GET_VALUE(4*163,x_af,y_af,x_a,y_a); b[163] = GET_VALUE(163*4 + 2,x_bf,y_bf,x_b,y_b);
            a[164] = GET_VALUE(4*164,x_af,y_af,x_a,y_a); b[164] = GET_VALUE(164*4 + 2,x_bf,y_bf,x_b,y_b);
            a[165] = GET_VALUE(4*165,x_af,y_af,x_a,y_a); b[165] = GET_VALUE(165*4 + 2,x_bf,y_bf,x_b,y_b);
            a[166] = GET_VALUE(4*166,x_af,y_af,x_a,y_a); b[166] = GET_VALUE(166*4 + 2,x_bf,y_bf,x_b,y_b);
            a[167] = GET_VALUE(4*167,x_af,y_af,x_a,y_a); b[167] = GET_VALUE(167*4 + 2,x_bf,y_bf,x_b,y_b);
            a[168] = GET_VALUE(4*168,x_af,y_af,x_a,y_a); b[168] = GET_VALUE(168*4 + 2,x_bf,y_bf,x_b,y_b);
            a[169] = GET_VALUE(4*169,x_af,y_af,x_a,y_a); b[169] = GET_VALUE(169*4 + 2,x_bf,y_bf,x_b,y_b);
            a[170] = GET_VALUE(4*170,x_af,y_af,x_a,y_a); b[170] = GET_VALUE(170*4 + 2,x_bf,y_bf,x_b,y_b);
            a[171] = GET_VALUE(4*171,x_af,y_af,x_a,y_a); b[171] = GET_VALUE(171*4 + 2,x_bf,y_bf,x_b,y_b);
            a[172] = GET_VALUE(4*172,x_af,y_af,x_a,y_a); b[172] = GET_VALUE(172*4 + 2,x_bf,y_bf,x_b,y_b);
            a[173] = GET_VALUE(4*173,x_af,y_af,x_a,y_a); b[173] = GET_VALUE(173*4 + 2,x_bf,y_bf,x_b,y_b);
            a[174] = GET_VALUE(4*174,x_af,y_af,x_a,y_a); b[174] = GET_VALUE(174*4 + 2,x_bf,y_bf,x_b,y_b);
            a[175] = GET_VALUE(4*175,x_af,y_af,x_a,y_a); b[175] = GET_VALUE(175*4 + 2,x_bf,y_bf,x_b,y_b);
            a[176] = GET_VALUE(4*176,x_af,y_af,x_a,y_a); b[176] = GET_VALUE(176*4 + 2,x_bf,y_bf,x_b,y_b);
            a[177] = GET_VALUE(4*177,x_af,y_af,x_a,y_a); b[177] = GET_VALUE(177*4 + 2,x_bf,y_bf,x_b,y_b);
            a[178] = GET_VALUE(4*178,x_af,y_af,x_a,y_a); b[178] = GET_VALUE(178*4 + 2,x_bf,y_bf,x_b,y_b);
            a[179] = GET_VALUE(4*179,x_af,y_af,x_a,y_a); b[179] = GET_VALUE(179*4 + 2,x_bf,y_bf,x_b,y_b);
            a[180] = GET_VALUE(4*180,x_af,y_af,x_a,y_a); b[180] = GET_VALUE(180*4 + 2,x_bf,y_bf,x_b,y_b);
            a[181] = GET_VALUE(4*181,x_af,y_af,x_a,y_a); b[181] = GET_VALUE(181*4 + 2,x_bf,y_bf,x_b,y_b);
            a[182] = GET_VALUE(4*182,x_af,y_af,x_a,y_a); b[182] = GET_VALUE(182*4 + 2,x_bf,y_bf,x_b,y_b);
            a[183] = GET_VALUE(4*183,x_af,y_af,x_a,y_a); b[183] = GET_VALUE(183*4 + 2,x_bf,y_bf,x_b,y_b);
            a[184] = GET_VALUE(4*184,x_af,y_af,x_a,y_a); b[184] = GET_VALUE(184*4 + 2,x_bf,y_bf,x_b,y_b);
            a[185] = GET_VALUE(4*185,x_af,y_af,x_a,y_a); b[185] = GET_VALUE(185*4 + 2,x_bf,y_bf,x_b,y_b);
            a[186] = GET_VALUE(4*186,x_af,y_af,x_a,y_a); b[186] = GET_VALUE(186*4 + 2,x_bf,y_bf,x_b,y_b);
            a[187] = GET_VALUE(4*187,x_af,y_af,x_a,y_a); b[187] = GET_VALUE(187*4 + 2,x_bf,y_bf,x_b,y_b);
            a[188] = GET_VALUE(4*188,x_af,y_af,x_a,y_a); b[188] = GET_VALUE(188*4 + 2,x_bf,y_bf,x_b,y_b);
            a[189] = GET_VALUE(4*189,x_af,y_af,x_a,y_a); b[189] = GET_VALUE(189*4 + 2,x_bf,y_bf,x_b,y_b);
            a[190] = GET_VALUE(4*190,x_af,y_af,x_a,y_a); b[190] = GET_VALUE(190*4 + 2,x_bf,y_bf,x_b,y_b);
            a[191] = GET_VALUE(4*191,x_af,y_af,x_a,y_a); b[191] = GET_VALUE(191*4 + 2,x_bf,y_bf,x_b,y_b);
            a[192] = GET_VALUE(4*192,x_af,y_af,x_a,y_a); b[192] = GET_VALUE(192*4 + 2,x_bf,y_bf,x_b,y_b);
            a[193] = GET_VALUE(4*193,x_af,y_af,x_a,y_a); b[193] = GET_VALUE(193*4 + 2,x_bf,y_bf,x_b,y_b);
            a[194] = GET_VALUE(4*194,x_af,y_af,x_a,y_a); b[194] = GET_VALUE(194*4 + 2,x_bf,y_bf,x_b,y_b);
            a[195] = GET_VALUE(4*195,x_af,y_af,x_a,y_a); b[195] = GET_VALUE(195*4 + 2,x_bf,y_bf,x_b,y_b);
            a[196] = GET_VALUE(4*196,x_af,y_af,x_a,y_a); b[196] = GET_VALUE(196*4 + 2,x_bf,y_bf,x_b,y_b);
            a[197] = GET_VALUE(4*197,x_af,y_af,x_a,y_a); b[197] = GET_VALUE(197*4 + 2,x_bf,y_bf,x_b,y_b);
            a[198] = GET_VALUE(4*198,x_af,y_af,x_a,y_a); b[198] = GET_VALUE(198*4 + 2,x_bf,y_bf,x_b,y_b);
            a[199] = GET_VALUE(4*199,x_af,y_af,x_a,y_a); b[199] = GET_VALUE(199*4 + 2,x_bf,y_bf,x_b,y_b);
            a[200] = GET_VALUE(4*200,x_af,y_af,x_a,y_a); b[200] = GET_VALUE(200*4 + 2,x_bf,y_bf,x_b,y_b);
            a[201] = GET_VALUE(4*201,x_af,y_af,x_a,y_a); b[201] = GET_VALUE(201*4 + 2,x_bf,y_bf,x_b,y_b);
            a[202] = GET_VALUE(4*202,x_af,y_af,x_a,y_a); b[202] = GET_VALUE(202*4 + 2,x_bf,y_bf,x_b,y_b);
            a[203] = GET_VALUE(4*203,x_af,y_af,x_a,y_a); b[203] = GET_VALUE(203*4 + 2,x_bf,y_bf,x_b,y_b);
            a[204] = GET_VALUE(4*204,x_af,y_af,x_a,y_a); b[204] = GET_VALUE(204*4 + 2,x_bf,y_bf,x_b,y_b);
            a[205] = GET_VALUE(4*205,x_af,y_af,x_a,y_a); b[205] = GET_VALUE(205*4 + 2,x_bf,y_bf,x_b,y_b);
            a[206] = GET_VALUE(4*206,x_af,y_af,x_a,y_a); b[206] = GET_VALUE(206*4 + 2,x_bf,y_bf,x_b,y_b);
            a[207] = GET_VALUE(4*207,x_af,y_af,x_a,y_a); b[207] = GET_VALUE(207*4 + 2,x_bf,y_bf,x_b,y_b);
            a[208] = GET_VALUE(4*208,x_af,y_af,x_a,y_a); b[208] = GET_VALUE(208*4 + 2,x_bf,y_bf,x_b,y_b);
            a[209] = GET_VALUE(4*209,x_af,y_af,x_a,y_a); b[209] = GET_VALUE(209*4 + 2,x_bf,y_bf,x_b,y_b);
            a[210] = GET_VALUE(4*210,x_af,y_af,x_a,y_a); b[210] = GET_VALUE(210*4 + 2,x_bf,y_bf,x_b,y_b);
            a[211] = GET_VALUE(4*211,x_af,y_af,x_a,y_a); b[211] = GET_VALUE(211*4 + 2,x_bf,y_bf,x_b,y_b);
            a[212] = GET_VALUE(4*212,x_af,y_af,x_a,y_a); b[212] = GET_VALUE(212*4 + 2,x_bf,y_bf,x_b,y_b);
            a[213] = GET_VALUE(4*213,x_af,y_af,x_a,y_a); b[213] = GET_VALUE(213*4 + 2,x_bf,y_bf,x_b,y_b);
            a[214] = GET_VALUE(4*214,x_af,y_af,x_a,y_a); b[214] = GET_VALUE(214*4 + 2,x_bf,y_bf,x_b,y_b);
            a[215] = GET_VALUE(4*215,x_af,y_af,x_a,y_a); b[215] = GET_VALUE(215*4 + 2,x_bf,y_bf,x_b,y_b);
            a[216] = GET_VALUE(4*216,x_af,y_af,x_a,y_a); b[216] = GET_VALUE(216*4 + 2,x_bf,y_bf,x_b,y_b);
            a[217] = GET_VALUE(4*217,x_af,y_af,x_a,y_a); b[217] = GET_VALUE(217*4 + 2,x_bf,y_bf,x_b,y_b);
            a[218] = GET_VALUE(4*218,x_af,y_af,x_a,y_a); b[218] = GET_VALUE(218*4 + 2,x_bf,y_bf,x_b,y_b);
            a[219] = GET_VALUE(4*219,x_af,y_af,x_a,y_a); b[219] = GET_VALUE(219*4 + 2,x_bf,y_bf,x_b,y_b);
            a[220] = GET_VALUE(4*220,x_af,y_af,x_a,y_a); b[220] = GET_VALUE(220*4 + 2,x_bf,y_bf,x_b,y_b);
            a[221] = GET_VALUE(4*221,x_af,y_af,x_a,y_a); b[221] = GET_VALUE(221*4 + 2,x_bf,y_bf,x_b,y_b);
            a[222] = GET_VALUE(4*222,x_af,y_af,x_a,y_a); b[222] = GET_VALUE(222*4 + 2,x_bf,y_bf,x_b,y_b);
            a[223] = GET_VALUE(4*223,x_af,y_af,x_a,y_a); b[223] = GET_VALUE(223*4 + 2,x_bf,y_bf,x_b,y_b);
            a[224] = GET_VALUE(4*224,x_af,y_af,x_a,y_a); b[224] = GET_VALUE(224*4 + 2,x_bf,y_bf,x_b,y_b);
            a[225] = GET_VALUE(4*225,x_af,y_af,x_a,y_a); b[225] = GET_VALUE(225*4 + 2,x_bf,y_bf,x_b,y_b);
            a[226] = GET_VALUE(4*226,x_af,y_af,x_a,y_a); b[226] = GET_VALUE(226*4 + 2,x_bf,y_bf,x_b,y_b);
            a[227] = GET_VALUE(4*227,x_af,y_af,x_a,y_a); b[227] = GET_VALUE(227*4 + 2,x_bf,y_bf,x_b,y_b);
            a[228] = GET_VALUE(4*228,x_af,y_af,x_a,y_a); b[228] = GET_VALUE(228*4 + 2,x_bf,y_bf,x_b,y_b);
            a[229] = GET_VALUE(4*229,x_af,y_af,x_a,y_a); b[229] = GET_VALUE(229*4 + 2,x_bf,y_bf,x_b,y_b);
            a[230] = GET_VALUE(4*230,x_af,y_af,x_a,y_a); b[230] = GET_VALUE(230*4 + 2,x_bf,y_bf,x_b,y_b);
            a[231] = GET_VALUE(4*231,x_af,y_af,x_a,y_a); b[231] = GET_VALUE(231*4 + 2,x_bf,y_bf,x_b,y_b);
            a[232] = GET_VALUE(4*232,x_af,y_af,x_a,y_a); b[232] = GET_VALUE(232*4 + 2,x_bf,y_bf,x_b,y_b);
            a[233] = GET_VALUE(4*233,x_af,y_af,x_a,y_a); b[233] = GET_VALUE(233*4 + 2,x_bf,y_bf,x_b,y_b);
            a[234] = GET_VALUE(4*234,x_af,y_af,x_a,y_a); b[234] = GET_VALUE(234*4 + 2,x_bf,y_bf,x_b,y_b);
            a[235] = GET_VALUE(4*235,x_af,y_af,x_a,y_a); b[235] = GET_VALUE(235*4 + 2,x_bf,y_bf,x_b,y_b);
            a[236] = GET_VALUE(4*236,x_af,y_af,x_a,y_a); b[236] = GET_VALUE(236*4 + 2,x_bf,y_bf,x_b,y_b);
            a[237] = GET_VALUE(4*237,x_af,y_af,x_a,y_a); b[237] = GET_VALUE(237*4 + 2,x_bf,y_bf,x_b,y_b);
            a[238] = GET_VALUE(4*238,x_af,y_af,x_a,y_a); b[238] = GET_VALUE(238*4 + 2,x_bf,y_bf,x_b,y_b);
            a[239] = GET_VALUE(4*239,x_af,y_af,x_a,y_a); b[239] = GET_VALUE(239*4 + 2,x_bf,y_bf,x_b,y_b);
            a[240] = GET_VALUE(4*240,x_af,y_af,x_a,y_a); b[240] = GET_VALUE(240*4 + 2,x_bf,y_bf,x_b,y_b);
            a[241] = GET_VALUE(4*241,x_af,y_af,x_a,y_a); b[241] = GET_VALUE(241*4 + 2,x_bf,y_bf,x_b,y_b);
            a[242] = GET_VALUE(4*242,x_af,y_af,x_a,y_a); b[242] = GET_VALUE(242*4 + 2,x_bf,y_bf,x_b,y_b);
            a[243] = GET_VALUE(4*243,x_af,y_af,x_a,y_a); b[243] = GET_VALUE(243*4 + 2,x_bf,y_bf,x_b,y_b);
            a[244] = GET_VALUE(4*244,x_af,y_af,x_a,y_a); b[244] = GET_VALUE(244*4 + 2,x_bf,y_bf,x_b,y_b);
            a[245] = GET_VALUE(4*245,x_af,y_af,x_a,y_a); b[245] = GET_VALUE(245*4 + 2,x_bf,y_bf,x_b,y_b);
            a[246] = GET_VALUE(4*246,x_af,y_af,x_a,y_a); b[246] = GET_VALUE(246*4 + 2,x_bf,y_bf,x_b,y_b);
            a[247] = GET_VALUE(4*247,x_af,y_af,x_a,y_a); b[247] = GET_VALUE(247*4 + 2,x_bf,y_bf,x_b,y_b);
            a[248] = GET_VALUE(4*248,x_af,y_af,x_a,y_a); b[248] = GET_VALUE(248*4 + 2,x_bf,y_bf,x_b,y_b);
            a[249] = GET_VALUE(4*249,x_af,y_af,x_a,y_a); b[249] = GET_VALUE(249*4 + 2,x_bf,y_bf,x_b,y_b);
            a[250] = GET_VALUE(4*250,x_af,y_af,x_a,y_a); b[250] = GET_VALUE(250*4 + 2,x_bf,y_bf,x_b,y_b);
            a[251] = GET_VALUE(4*251,x_af,y_af,x_a,y_a); b[251] = GET_VALUE(251*4 + 2,x_bf,y_bf,x_b,y_b);
            a[252] = GET_VALUE(4*252,x_af,y_af,x_a,y_a); b[252] = GET_VALUE(252*4 + 2,x_bf,y_bf,x_b,y_b);
            a[253] = GET_VALUE(4*253,x_af,y_af,x_a,y_a); b[253] = GET_VALUE(253*4 + 2,x_bf,y_bf,x_b,y_b);
            a[254] = GET_VALUE(4*254,x_af,y_af,x_a,y_a); b[254] = GET_VALUE(254*4 + 2,x_bf,y_bf,x_b,y_b);
            a[255] = GET_VALUE(4*255,x_af,y_af,x_a,y_a); b[255] = GET_VALUE(255*4 + 2,x_bf,y_bf,x_b,y_b);
            
            int32_t f[16];
            f[0] =_mm_movemask_epi8(_mm_cmpgt_epi8(_mm_load_si128((__m128i const*)(a+16*0)),_mm_load_si128((__m128i const*) (b+16*0))));
            f[1] =_mm_movemask_epi8(_mm_cmpgt_epi8(_mm_load_si128((__m128i const*)(a+16*1)),_mm_load_si128((__m128i const*) (b+16*1))));
            f[2] =_mm_movemask_epi8(_mm_cmpgt_epi8(_mm_load_si128((__m128i const*)(a+16*2)),_mm_load_si128((__m128i const*) (b+16*2))));
            f[3] =_mm_movemask_epi8(_mm_cmpgt_epi8(_mm_load_si128((__m128i const*)(a+16*3)),_mm_load_si128((__m128i const*) (b+16*3))));
            f[4] =_mm_movemask_epi8(_mm_cmpgt_epi8(_mm_load_si128((__m128i const*)(a+16*4)),_mm_load_si128((__m128i const*) (b+16*4))));
            f[5] =_mm_movemask_epi8(_mm_cmpgt_epi8(_mm_load_si128((__m128i const*)(a+16*5)),_mm_load_si128((__m128i const*) (b+16*5))));
            f[6] =_mm_movemask_epi8(_mm_cmpgt_epi8(_mm_load_si128((__m128i const*)(a+16*6)),_mm_load_si128((__m128i const*) (b+16*6))));
            f[7] =_mm_movemask_epi8(_mm_cmpgt_epi8(_mm_load_si128((__m128i const*)(a+16*7)),_mm_load_si128((__m128i const*) (b+16*7))));
            f[8] =_mm_movemask_epi8(_mm_cmpgt_epi8(_mm_load_si128((__m128i const*)(a+16*8)),_mm_load_si128((__m128i const*) (b+16*8))));
            f[9] =_mm_movemask_epi8(_mm_cmpgt_epi8(_mm_load_si128((__m128i const*)(a+16*9)),_mm_load_si128((__m128i const*) (b+16*9))));
            f[10] =_mm_movemask_epi8(_mm_cmpgt_epi8(_mm_load_si128((__m128i const*)(a+16*10)),_mm_load_si128((__m128i const*) (b+16*10))));
            f[11] =_mm_movemask_epi8(_mm_cmpgt_epi8(_mm_load_si128((__m128i const*)(a+16*11)),_mm_load_si128((__m128i const*) (b+16*11))));
            f[12] =_mm_movemask_epi8(_mm_cmpgt_epi8(_mm_load_si128((__m128i const*)(a+16*12)),_mm_load_si128((__m128i const*) (b+16*12))));
            f[13] =_mm_movemask_epi8(_mm_cmpgt_epi8(_mm_load_si128((__m128i const*)(a+16*13)),_mm_load_si128((__m128i const*) (b+16*13))));
            f[14] =_mm_movemask_epi8(_mm_cmpgt_epi8(_mm_load_si128((__m128i const*)(a+16*14)),_mm_load_si128((__m128i const*) (b+16*14))));
            f[15] =_mm_movemask_epi8(_mm_cmpgt_epi8(_mm_load_si128((__m128i const*)(a+16*15)),_mm_load_si128((__m128i const*) (b+16*15))));
            
            //int32_t g[16] __attribute__((aligned(32)));

            bd[j*n_rows_bd] = ((int64_t)(f[0])) + ((int64_t)(f[1]) << 1) + ((int64_t)(f[2])<<2) + ((int64_t)(f[3])<<3);
        }
    }
}


int BRIEF::diag_length_pattern = 17;
int BRIEF::gaussian_bit_pattern_31[256 * 4] = {
        8, -3, 9, 5,/*mean (0), correlation (0)*/
        4, 2, 7, -12,/*mean (1.12461e-05), correlation (0.0437584)*/
        -11, 9, -8, 2,/*mean (3.37382e-05), correlation (0.0617409)*/
        7, -12, 12, -13,/*mean (5.62303e-05), correlation (0.0636977)*/
        2, -13, 2, 12,/*mean (0.000134953), correlation (0.085099)*/
        1, -7, 1, 6,/*mean (0.000528565), correlation (0.0857175)*/
        -2, -10, -2, -4,/*mean (0.0188821), correlation (0.0985774)*/
        -13, -13, -11, -8,/*mean (0.0363135), correlation (0.0899616)*/
        -13, -3, -12, -9,/*mean (0.121806), correlation (0.099849)*/
        10, 4, 11, 9,/*mean (0.122065), correlation (0.093285)*/
        -13, -8, -8, -9,/*mean (0.162787), correlation (0.0942748)*/
        -11, 7, -9, 12,/*mean (0.21561), correlation (0.0974438)*/
        7, 7, 12, 6,/*mean (0.160583), correlation (0.130064)*/
        -4, -5, -3, 0,/*mean (0.228171), correlation (0.132998)*/
        -13, 2, -12, -3,/*mean (0.00997526), correlation (0.145926)*/
        -9, 0, -7, 5,/*mean (0.198234), correlation (0.143636)*/
        12, -6, 12, -1,/*mean (0.0676226), correlation (0.16689)*/
        -3, 6, -2, 12,/*mean (0.166847), correlation (0.171682)*/
        -6, -13, -4, -8,/*mean (0.101215), correlation (0.179716)*/
        11, -13, 12, -8,/*mean (0.200641), correlation (0.192279)*/
        4, 7, 5, 1,/*mean (0.205106), correlation (0.186848)*/
        5, -3, 10, -3,/*mean (0.234908), correlation (0.192319)*/
        3, -7, 6, 12,/*mean (0.0709964), correlation (0.210872)*/
        -8, -7, -6, -2,/*mean (0.0939834), correlation (0.212589)*/
        -2, 11, -1, -10,/*mean (0.127778), correlation (0.20866)*/
        -13, 12, -8, 10,/*mean (0.14783), correlation (0.206356)*/
        -7, 3, -5, -3,/*mean (0.182141), correlation (0.198942)*/
        -4, 2, -3, 7,/*mean (0.188237), correlation (0.21384)*/
        -10, -12, -6, 11,/*mean (0.14865), correlation (0.23571)*/
        5, -12, 6, -7,/*mean (0.222312), correlation (0.23324)*/
        5, -6, 7, -1,/*mean (0.229082), correlation (0.23389)*/
        1, 0, 4, -5,/*mean (0.241577), correlation (0.215286)*/
        9, 11, 11, -13,/*mean (0.00338507), correlation (0.251373)*/
        4, 7, 4, 12,/*mean (0.131005), correlation (0.257622)*/
        2, -1, 4, 4,/*mean (0.152755), correlation (0.255205)*/
        -4, -12, -2, 7,/*mean (0.182771), correlation (0.244867)*/
        -8, -5, -7, -10,/*mean (0.186898), correlation (0.23901)*/
        4, 11, 9, 12,/*mean (0.226226), correlation (0.258255)*/
        0, -8, 1, -13,/*mean (0.0897886), correlation (0.274827)*/
        -13, -2, -8, 2,/*mean (0.148774), correlation (0.28065)*/
        -3, -2, -2, 3,/*mean (0.153048), correlation (0.283063)*/
        -6, 9, -4, -9,/*mean (0.169523), correlation (0.278248)*/
        8, 12, 10, 7,/*mean (0.225337), correlation (0.282851)*/
        0, 9, 1, 3,/*mean (0.226687), correlation (0.278734)*/
        7, -5, 11, -10,/*mean (0.00693882), correlation (0.305161)*/
        -13, -6, -11, 0,/*mean (0.0227283), correlation (0.300181)*/
        10, 7, 12, 1,/*mean (0.125517), correlation (0.31089)*/
        -6, -3, -6, 12,/*mean (0.131748), correlation (0.312779)*/
        10, -9, 12, -4,/*mean (0.144827), correlation (0.292797)*/
        -13, 8, -8, -12,/*mean (0.149202), correlation (0.308918)*/
        -13, 0, -8, -4,/*mean (0.160909), correlation (0.310013)*/
        3, 3, 7, 8,/*mean (0.177755), correlation (0.309394)*/
        5, 7, 10, -7,/*mean (0.212337), correlation (0.310315)*/
        -1, 7, 1, -12,/*mean (0.214429), correlation (0.311933)*/
        3, -10, 5, 6,/*mean (0.235807), correlation (0.313104)*/
        2, -4, 3, -10,/*mean (0.00494827), correlation (0.344948)*/
        -13, 0, -13, 5,/*mean (0.0549145), correlation (0.344675)*/
        -13, -7, -12, 12,/*mean (0.103385), correlation (0.342715)*/
        -13, 3, -11, 8,/*mean (0.134222), correlation (0.322922)*/
        -7, 12, -4, 7,/*mean (0.153284), correlation (0.337061)*/
        6, -10, 12, 8,/*mean (0.154881), correlation (0.329257)*/
        -9, -1, -7, -6,/*mean (0.200967), correlation (0.33312)*/
        -2, -5, 0, 12,/*mean (0.201518), correlation (0.340635)*/
        -12, 5, -7, 5,/*mean (0.207805), correlation (0.335631)*/
        3, -10, 8, -13,/*mean (0.224438), correlation (0.34504)*/
        -7, -7, -4, 5,/*mean (0.239361), correlation (0.338053)*/
        -3, -2, -1, -7,/*mean (0.240744), correlation (0.344322)*/
        2, 9, 5, -11,/*mean (0.242949), correlation (0.34145)*/
        -11, -13, -5, -13,/*mean (0.244028), correlation (0.336861)*/
        -1, 6, 0, -1,/*mean (0.247571), correlation (0.343684)*/
        5, -3, 5, 2,/*mean (0.000697256), correlation (0.357265)*/
        -4, -13, -4, 12,/*mean (0.00213675), correlation (0.373827)*/
        -9, -6, -9, 6,/*mean (0.0126856), correlation (0.373938)*/
        -12, -10, -8, -4,/*mean (0.0152497), correlation (0.364237)*/
        10, 2, 12, -3,/*mean (0.0299933), correlation (0.345292)*/
        7, 12, 12, 12,/*mean (0.0307242), correlation (0.366299)*/
        -7, -13, -6, 5,/*mean (0.0534975), correlation (0.368357)*/
        -4, 9, -3, 4,/*mean (0.099865), correlation (0.372276)*/
        7, -1, 12, 2,/*mean (0.117083), correlation (0.364529)*/
        -7, 6, -5, 1,/*mean (0.126125), correlation (0.369606)*/
        -13, 11, -12, 5,/*mean (0.130364), correlation (0.358502)*/
        -3, 7, -2, -6,/*mean (0.131691), correlation (0.375531)*/
        7, -8, 12, -7,/*mean (0.160166), correlation (0.379508)*/
        -13, -7, -11, -12,/*mean (0.167848), correlation (0.353343)*/
        1, -3, 12, 12,/*mean (0.183378), correlation (0.371916)*/
        2, -6, 3, 0,/*mean (0.228711), correlation (0.371761)*/
        -4, 3, -2, -13,/*mean (0.247211), correlation (0.364063)*/
        -1, -13, 1, 9,/*mean (0.249325), correlation (0.378139)*/
        7, 1, 8, -6,/*mean (0.000652272), correlation (0.411682)*/
        1, -1, 3, 12,/*mean (0.00248538), correlation (0.392988)*/
        9, 1, 12, 6,/*mean (0.0206815), correlation (0.386106)*/
        -1, -9, -1, 3,/*mean (0.0364485), correlation (0.410752)*/
        -13, -13, -10, 5,/*mean (0.0376068), correlation (0.398374)*/
        7, 7, 10, 12,/*mean (0.0424202), correlation (0.405663)*/
        12, -5, 12, 9,/*mean (0.0942645), correlation (0.410422)*/
        6, 3, 7, 11,/*mean (0.1074), correlation (0.413224)*/
        5, -13, 6, 10,/*mean (0.109256), correlation (0.408646)*/
        2, -12, 2, 3,/*mean (0.131691), correlation (0.416076)*/
        3, 8, 4, -6,/*mean (0.165081), correlation (0.417569)*/
        2, 6, 12, -13,/*mean (0.171874), correlation (0.408471)*/
        9, -12, 10, 3,/*mean (0.175146), correlation (0.41296)*/
        -8, 4, -7, 9,/*mean (0.183682), correlation (0.402956)*/
        -11, 12, -4, -6,/*mean (0.184672), correlation (0.416125)*/
        1, 12, 2, -8,/*mean (0.191487), correlation (0.386696)*/
        6, -9, 7, -4,/*mean (0.192668), correlation (0.394771)*/
        2, 3, 3, -2,/*mean (0.200157), correlation (0.408303)*/
        6, 3, 11, 0,/*mean (0.204588), correlation (0.411762)*/
        3, -3, 8, -8,/*mean (0.205904), correlation (0.416294)*/
        7, 8, 9, 3,/*mean (0.213237), correlation (0.409306)*/
        -11, -5, -6, -4,/*mean (0.243444), correlation (0.395069)*/
        -10, 11, -5, 10,/*mean (0.247672), correlation (0.413392)*/
        -5, -8, -3, 12,/*mean (0.24774), correlation (0.411416)*/
        -10, 5, -9, 0,/*mean (0.00213675), correlation (0.454003)*/
        8, -1, 12, -6,/*mean (0.0293635), correlation (0.455368)*/
        4, -6, 6, -11,/*mean (0.0404971), correlation (0.457393)*/
        -10, 12, -8, 7,/*mean (0.0481107), correlation (0.448364)*/
        4, -2, 6, 7,/*mean (0.050641), correlation (0.455019)*/
        -2, 0, -2, 12,/*mean (0.0525978), correlation (0.44338)*/
        -5, -8, -5, 2,/*mean (0.0629667), correlation (0.457096)*/
        7, -6, 10, 12,/*mean (0.0653846), correlation (0.445623)*/
        -9, -13, -8, -8,/*mean (0.0858749), correlation (0.449789)*/
        -5, -13, -5, -2,/*mean (0.122402), correlation (0.450201)*/
        8, -8, 9, -13,/*mean (0.125416), correlation (0.453224)*/
        -9, -11, -9, 0,/*mean (0.130128), correlation (0.458724)*/
        1, -8, 1, -2,/*mean (0.132467), correlation (0.440133)*/
        7, -4, 9, 1,/*mean (0.132692), correlation (0.454)*/
        -2, 1, -1, -4,/*mean (0.135695), correlation (0.455739)*/
        11, -6, 12, -11,/*mean (0.142904), correlation (0.446114)*/
        -12, -9, -6, 4,/*mean (0.146165), correlation (0.451473)*/
        3, 7, 7, 12,/*mean (0.147627), correlation (0.456643)*/
        5, 5, 10, 8,/*mean (0.152901), correlation (0.455036)*/
        0, -4, 2, 8,/*mean (0.167083), correlation (0.459315)*/
        -9, 12, -5, -13,/*mean (0.173234), correlation (0.454706)*/
        0, 7, 2, 12,/*mean (0.18312), correlation (0.433855)*/
        -1, 2, 1, 7,/*mean (0.185504), correlation (0.443838)*/
        5, 11, 7, -9,/*mean (0.185706), correlation (0.451123)*/
        3, 5, 6, -8,/*mean (0.188968), correlation (0.455808)*/
        -13, -4, -8, 9,/*mean (0.191667), correlation (0.459128)*/
        -5, 9, -3, -3,/*mean (0.193196), correlation (0.458364)*/
        -4, -7, -3, -12,/*mean (0.196536), correlation (0.455782)*/
        6, 5, 8, 0,/*mean (0.1972), correlation (0.450481)*/
        -7, 6, -6, 12,/*mean (0.199438), correlation (0.458156)*/
        -13, 6, -5, -2,/*mean (0.211224), correlation (0.449548)*/
        1, -10, 3, 10,/*mean (0.211718), correlation (0.440606)*/
        4, 1, 8, -4,/*mean (0.213034), correlation (0.443177)*/
        -2, -2, 2, -13,/*mean (0.234334), correlation (0.455304)*/
        2, -12, 12, 12,/*mean (0.235684), correlation (0.443436)*/
        -2, -13, 0, -6,/*mean (0.237674), correlation (0.452525)*/
        4, 1, 9, 3,/*mean (0.23962), correlation (0.444824)*/
        -6, -10, -3, -5,/*mean (0.248459), correlation (0.439621)*/
        -3, -13, -1, 1,/*mean (0.249505), correlation (0.456666)*/
        7, 5, 12, -11,/*mean (0.00119208), correlation (0.495466)*/
        4, -2, 5, -7,/*mean (0.00372245), correlation (0.484214)*/
        -13, 9, -9, -5,/*mean (0.00741116), correlation (0.499854)*/
        7, 1, 8, 6,/*mean (0.0208952), correlation (0.499773)*/
        7, -8, 7, 6,/*mean (0.0220085), correlation (0.501609)*/
        -7, -4, -7, 1,/*mean (0.0233806), correlation (0.496568)*/
        -8, 11, -7, -8,/*mean (0.0236505), correlation (0.489719)*/
        -13, 6, -12, -8,/*mean (0.0268781), correlation (0.503487)*/
        2, 4, 3, 9,/*mean (0.0323324), correlation (0.501938)*/
        10, -5, 12, 3,/*mean (0.0399235), correlation (0.494029)*/
        -6, -5, -6, 7,/*mean (0.0420153), correlation (0.486579)*/
        8, -3, 9, -8,/*mean (0.0548021), correlation (0.484237)*/
        2, -12, 2, 8,/*mean (0.0616622), correlation (0.496642)*/
        -11, -2, -10, 3,/*mean (0.0627755), correlation (0.498563)*/
        -12, -13, -7, -9,/*mean (0.0829622), correlation (0.495491)*/
        -11, 0, -10, -5,/*mean (0.0843342), correlation (0.487146)*/
        5, -3, 11, 8,/*mean (0.0929937), correlation (0.502315)*/
        -2, -13, -1, 12,/*mean (0.113327), correlation (0.48941)*/
        -1, -8, 0, 9,/*mean (0.132119), correlation (0.467268)*/
        -13, -11, -12, -5,/*mean (0.136269), correlation (0.498771)*/
        -10, -2, -10, 11,/*mean (0.142173), correlation (0.498714)*/
        -3, 9, -2, -13,/*mean (0.144141), correlation (0.491973)*/
        2, -3, 3, 2,/*mean (0.14892), correlation (0.500782)*/
        -9, -13, -4, 0,/*mean (0.150371), correlation (0.498211)*/
        -4, 6, -3, -10,/*mean (0.152159), correlation (0.495547)*/
        -4, 12, -2, -7,/*mean (0.156152), correlation (0.496925)*/
        -6, -11, -4, 9,/*mean (0.15749), correlation (0.499222)*/
        6, -3, 6, 11,/*mean (0.159211), correlation (0.503821)*/
        -13, 11, -5, 5,/*mean (0.162427), correlation (0.501907)*/
        11, 11, 12, 6,/*mean (0.16652), correlation (0.497632)*/
        7, -5, 12, -2,/*mean (0.169141), correlation (0.484474)*/
        -1, 12, 0, 7,/*mean (0.169456), correlation (0.495339)*/
        -4, -8, -3, -2,/*mean (0.171457), correlation (0.487251)*/
        -7, 1, -6, 7,/*mean (0.175), correlation (0.500024)*/
        -13, -12, -8, -13,/*mean (0.175866), correlation (0.497523)*/
        -7, -2, -6, -8,/*mean (0.178273), correlation (0.501854)*/
        -8, 5, -6, -9,/*mean (0.181107), correlation (0.494888)*/
        -5, -1, -4, 5,/*mean (0.190227), correlation (0.482557)*/
        -13, 7, -8, 10,/*mean (0.196739), correlation (0.496503)*/
        1, 5, 5, -13,/*mean (0.19973), correlation (0.499759)*/
        1, 0, 10, -13,/*mean (0.204465), correlation (0.49873)*/
        9, 12, 10, -1,/*mean (0.209334), correlation (0.49063)*/
        5, -8, 10, -9,/*mean (0.211134), correlation (0.503011)*/
        -1, 11, 1, -13,/*mean (0.212), correlation (0.499414)*/
        -9, -3, -6, 2,/*mean (0.212168), correlation (0.480739)*/
        -1, -10, 1, 12,/*mean (0.212731), correlation (0.502523)*/
        -13, 1, -8, -10,/*mean (0.21327), correlation (0.489786)*/
        8, -11, 10, -6,/*mean (0.214159), correlation (0.488246)*/
        2, -13, 3, -6,/*mean (0.216993), correlation (0.50287)*/
        7, -13, 12, -9,/*mean (0.223639), correlation (0.470502)*/
        -10, -10, -5, -7,/*mean (0.224089), correlation (0.500852)*/
        -10, -8, -8, -13,/*mean (0.228666), correlation (0.502629)*/
        4, -6, 8, 5,/*mean (0.22906), correlation (0.498305)*/
        3, 12, 8, -13,/*mean (0.233378), correlation (0.503825)*/
        -4, 2, -3, -3,/*mean (0.234323), correlation (0.476692)*/
        5, -13, 10, -12,/*mean (0.236392), correlation (0.475462)*/
        4, -13, 5, -1,/*mean (0.236842), correlation (0.504132)*/
        -9, 9, -4, 3,/*mean (0.236977), correlation (0.497739)*/
        0, 3, 3, -9,/*mean (0.24314), correlation (0.499398)*/
        -12, 1, -6, 1,/*mean (0.243297), correlation (0.489447)*/
        3, 2, 4, -8,/*mean (0.00155196), correlation (0.553496)*/
        -10, -10, -10, 9,/*mean (0.00239541), correlation (0.54297)*/
        8, -13, 12, 12,/*mean (0.0034413), correlation (0.544361)*/
        -8, -12, -6, -5,/*mean (0.003565), correlation (0.551225)*/
        2, 2, 3, 7,/*mean (0.00835583), correlation (0.55285)*/
        10, 6, 11, -8,/*mean (0.00885065), correlation (0.540913)*/
        6, 8, 8, -12,/*mean (0.0101552), correlation (0.551085)*/
        -7, 10, -6, 5,/*mean (0.0102227), correlation (0.533635)*/
        -3, -9, -3, 9,/*mean (0.0110211), correlation (0.543121)*/
        -1, -13, -1, 5,/*mean (0.0113473), correlation (0.550173)*/
        -3, -7, -3, 4,/*mean (0.0140913), correlation (0.554774)*/
        -8, -2, -8, 3,/*mean (0.017049), correlation (0.55461)*/
        4, 2, 12, 12,/*mean (0.01778), correlation (0.546921)*/
        2, -5, 3, 11,/*mean (0.0224022), correlation (0.549667)*/
        6, -9, 11, -13,/*mean (0.029161), correlation (0.546295)*/
        3, -1, 7, 12,/*mean (0.0303081), correlation (0.548599)*/
        11, -1, 12, 4,/*mean (0.0355151), correlation (0.523943)*/
        -3, 0, -3, 6,/*mean (0.0417904), correlation (0.543395)*/
        4, -11, 4, 12,/*mean (0.0487292), correlation (0.542818)*/
        2, -4, 2, 1,/*mean (0.0575124), correlation (0.554888)*/
        -10, -6, -8, 1,/*mean (0.0594242), correlation (0.544026)*/
        -13, 7, -11, 1,/*mean (0.0597391), correlation (0.550524)*/
        -13, 12, -11, -13,/*mean (0.0608974), correlation (0.55383)*/
        6, 0, 11, -13,/*mean (0.065126), correlation (0.552006)*/
        0, -1, 1, 4,/*mean (0.074224), correlation (0.546372)*/
        -13, 3, -9, -2,/*mean (0.0808592), correlation (0.554875)*/
        -9, 8, -6, -3,/*mean (0.0883378), correlation (0.551178)*/
        -13, -6, -8, -2,/*mean (0.0901035), correlation (0.548446)*/
        5, -9, 8, 10,/*mean (0.0949843), correlation (0.554694)*/
        2, 7, 3, -9,/*mean (0.0994152), correlation (0.550979)*/
        -1, -6, -1, -1,/*mean (0.10045), correlation (0.552714)*/
        9, 5, 11, -2,/*mean (0.100686), correlation (0.552594)*/
        11, -3, 12, -8,/*mean (0.101091), correlation (0.532394)*/
        3, 0, 3, 5,/*mean (0.101147), correlation (0.525576)*/
        -1, 4, 0, 10,/*mean (0.105263), correlation (0.531498)*/
        3, -6, 4, 5,/*mean (0.110785), correlation (0.540491)*/
        -13, 0, -10, 5,/*mean (0.112798), correlation (0.536582)*/
        5, 8, 12, 11,/*mean (0.114181), correlation (0.555793)*/
        8, 9, 9, -6,/*mean (0.117431), correlation (0.553763)*/
        7, -4, 8, -12,/*mean (0.118522), correlation (0.553452)*/
        -10, 4, -10, 9,/*mean (0.12094), correlation (0.554785)*/
        7, 3, 12, 4,/*mean (0.122582), correlation (0.555825)*/
        9, -7, 10, -2,/*mean (0.124978), correlation (0.549846)*/
        7, 0, 12, -2,/*mean (0.127002), correlation (0.537452)*/
        -1, -6, 0, -11,/*mean (0.127148), correlation (0.547401)*/
};
