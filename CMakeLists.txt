cmake_minimum_required(VERSION 3.6)
project(argon_simulation)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES basic_msim.cpp)

add_executable(argon_simulation ${SOURCE_FILES})