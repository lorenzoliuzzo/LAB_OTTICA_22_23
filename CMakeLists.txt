cmake_minimum_required(VERSION 3.20)

project(lab_fisica_moderna VERSION 2.0 LANGUAGES CXX)

set(CMAKE_CXX_FLAGS "-std=c++23 -O3 -g --pedantic -Wall -Wextra -Wno-literal-suffix -fconcepts -flto -ftree-loop-vectorize -I/usr/local/include/scipp")

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/lab_fisica_moderna)