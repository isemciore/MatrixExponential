cmake_minimum_required(VERSION 3.3)
project(task2folder)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp r8mat_expm1.cpp r8lib.cpp Matrix.h r8mat_expm1.h )
add_executable(task2folder ${SOURCE_FILES})