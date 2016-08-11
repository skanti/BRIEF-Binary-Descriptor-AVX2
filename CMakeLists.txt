cmake_minimum_required(VERSION 3.0)
project(BRIEF CXX)

set(CMAKE_CXX_COMPILER "/usr/local/bin/icpc")
set(OPTIMIZATION_FLAGS "-O3 -qopenmp -fomit-frame-pointer -march=native -ansi-alias")
set(USER_DEFINES "-DSIZE_BITS_HAMING=64 -DN_THREADS=4 -DN_DIM_BINARYDESCRIPTOR=256")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPTIMIZATION_FLAGS} -Wall -std=c++14 ${USER_DEFINES}")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")


set(SOURCE_FILES_CPP
        src/BRIEF.cpp)

set(SOURCE_FILES_H
        src/BRIEF.h
        src/Types.h)


set(OPENCV_DIR /usr/local/Cellar/opencv3/3.1.0_3)

include_directories(${OPENCV_DIR}/include ./src)
link_directories(${OPENCV_DIR}/lib)
add_library(BRIEF ${SOURCE_FILES_CPP} ${SOURCE_FILES_H})

add_executable(BRIEFExample example/main.cpp example/Timer.h)
target_link_libraries(BRIEFExample BRIEF opencv_core opencv_imgcodecs opencv_highgui)

add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/src/BRIEF.cpp
        PRE_BUILD
        COMMAND m4 -DSIZE_BITS_HAMING=64 -DN_THREADS=4 -DN_DIM_BINARYDESCRIPTOR=256 ${CMAKE_CURRENT_SOURCE_DIR}/src/BRIEF.m4 > $@
        MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/src/BRIEF.m4
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT ${USER_DEFINES})

set_source_files_properties(src/BRIEF.cpp PROPERTIES GENERATED true)