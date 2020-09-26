cmake_minimum_required(VERSION 3.0.0)
project(StarStack)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w -g3 -ggdb -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w -g3 -ggdb -std=c++11 -O3 -march=native -fopenmp")

add_subdirectory(image_utils)
include_directories(image_utils)
add_executable(starstack main.cpp $<TARGET_OBJECTS:image_utils> astro.cpp cfm.cpp tests.cpp)