cmake_minimum_required(VERSION 3.0)
project(BRIEF CXX)

set(CMAKE_CXX_COMPILER "/usr/local/bin/icpc")
set(OPTIMIZATION_FLAGS "-O3 -qopenmp -fomit-frame-pointer -march=native -ansi-alias")
set(USER_DEFINES "-DN_THREADS=4")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPTIMIZATION_FLAGS} -Wall -std=c++14 ${USER_DEFINES}")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")

set(SOURCE_FILES_O
        build/MathKernels.o)

set(SOURCE_FILES_CPP
        src/BRIEF.cpp)

set(SOURCE_FILES_H
        src/BRIEF.h
        src/Types.h)


set(OPENCV_DIR /usr/local/Cellar/opencv3/3.1.0_3)

include_directories(${OPENCV_DIR}/include ./src)
link_directories(${OPENCV_DIR}/lib)
add_library(BRIEF ${SOURCE_FILES_CPP} ${SOURCE_FILES_H} ${SOURCE_FILES_O})

add_executable(BRIEFExample example/main.cpp example/Timer.h)
target_link_libraries(BRIEFExample BRIEF opencv_core opencv_imgcodecs opencv_highgui)

add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/src/BRIEF.cpp
        PRE_BUILD
        COMMAND m4 -DSIZE_BITS_HAMING=64 -DN_DIM_BINARYDESCRIPTOR=256 ${CMAKE_CURRENT_SOURCE_DIR}/src/BRIEF.m4 > $@
        MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/src/BRIEF.m4
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT ${USER_DEFINES})

add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/build/MathKernels.o
        PRE_BUILD
        COMMAND ispc -O2 --arch=x86-64 --target=avx2 --math-lib=fast ${CMAKE_CURRENT_SOURCE_DIR}/src/MathKernels.ispc -o $@
        MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/src/MathKernels.ispc
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT "Compiling ISPC files.")


set_source_files_properties(src/BRIEF.cpp PROPERTIES GENERATED true)
set_source_files_properties(build/MathKernels.o PROPERTIES EXTERNAL_OBJECT true GENERATED true)