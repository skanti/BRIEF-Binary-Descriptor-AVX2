cmake_minimum_required(VERSION 3.0)
project(BRIEF CXX)

# -> required libraries
find_package(OpenCV 2.4 REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})
# <-

# -> compiler requirements
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED 1)
# <-

set(CMAKE_CXX_COMPILER "/usr/bin/c++")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-result -O3 -mavx2 -msse4 -Wall -std=c++14 -DNDEBUG")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")

set(SOURCE_FILES_CPP
        src/BRIEF.cpp)

set(SOURCE_FILES_H
        src/BRIEF.h
        src/Types.h)


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
add_library(BRIEF ${SOURCE_FILES_CPP} ${SOURCE_FILES_H})

#opencv_core opencv_imgcodecs opencv_highgui 
add_executable(BRIEFExample example/main.cpp example/Timer.h)
target_link_libraries(BRIEFExample BRIEF ${OpenCV_LIBRARIES})

add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/src/BRIEF.h
        PRE_BUILD
        COMMAND m4 -DSIZE_BITS_HAMING=64 -DN_DIM_BINARYDESCRIPTOR=256 ${CMAKE_CURRENT_SOURCE_DIR}/src/BRIEF.m4 > ${CMAKE_CURRENT_SOURCE_DIR}/src/BRIEF.h
        MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/src/BRIEF.m4
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT "Creating files from m4.")

set_source_files_properties(src/BRIEF.h PROPERTIES GENERATED true)