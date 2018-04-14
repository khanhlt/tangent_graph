cmake_minimum_required(VERSION 3.9)
project(tangent_graph)

set(CMAKE_CXX_STANDARD 11)

add_executable(tangent_graph main.cpp geomathhelper/geo_math_helper.cc geometry_library/geo_lib.cc graph/o_graph.cc)