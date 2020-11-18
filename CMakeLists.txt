cmake_minimum_required(VERSION 2.6)

project(Packing3D)

set(CMAKE_CXX_FLAGS "-std=c++11")

find_package(JNI REQUIRED)

add_library(packing3d STATIC main.cpp)

add_executable(test_packing3d main.cpp)
