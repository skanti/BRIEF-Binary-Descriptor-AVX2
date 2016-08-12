# BRIEF-Binary-Descriptor-AVX2
A fast BRIEF Implementation for 256-dimensional 64-bit Binary Descriptors

## Description
This is a C++ method that allows you to calculate the BRIEF binary feature descriptors.

## Performance
This method is approx. 8x faster than OpenCV's ```BRIEF``` method. The reason is that I am using a 64-bit strategy
as opposed to the 8-bit strategy of OpenCV. Furthermore, my code is aggressively optimized. That is, the kernel contains
mostly packed avx2 instructions.

#### Timing for a set of 500 features for 256-dimensional binary descriptors:
- OpenCV: 2.5 ms
- this version: 0.31 ms

## Features
- ```32-byte aligned plain arrays``` instead ```stl vectors```
- Gaussian pattern is divided into 4 256x1 contiguous arrays instead of a 256x4 matrices. Allows vectorization.
- Complete unrolling of the ```i = 0 ... 256``` loop
- AVX2 masking (```_mm256_movemask_epi8```) and comparing (```_mm256_cmpeq_epi8```) and load-n-store removes need of bitshifting

## Photos
#### PHOTO1: Assembly output shows nice and consistent packed instructions
<a href="http://tinypic.com?ref=34z0gh1" target="_blank"><img src="http://i65.tinypic.com/34z0gh1.jpg" border="0" ></a>
#### PHOTO2: Intel intrinsics handle gathering and ordering of bits in a packed way:
<a href="http://tinypic.com?ref=im2xkp" target="_blank"><img src="http://i63.tinypic.com/im2xkp.jpg" border="0"></a>

## To-Do
- [x] Unroll loops by pre-processor
- [ ] Try to speed up memory-bound bottleneck by AVX2 gather/scatter instructions

## Dependencies
- m4 (macro pre-processor)
- Intel ispc SIMD compiler

## How-To Install
Have a look at the CMake file. Just run it with your individual paths.

## Contact
Feel free to contact me if you have questions or just want to chat about it.
