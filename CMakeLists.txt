cmake_minimum_required(VERSION 2.6)

project(Packing3D)

set(CMAKE_CXX_FLAGS "-std=c++11")

find_package(JNI REQUIRED)

include_directories(${JNI_INCLUDE_DIRS})
# link_directories(deps)
link_libraries(${JNI_LIBRARIES})

add_library(packing3d STATIC main.cpp)

add_executable(test_packing3d main.cpp)
