# 设置项目名称
cmake_minimum_required(VERSION 3.10)
project(Graph)

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 14)

# 设置头文件目录
include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(Graph ${PROJECT_SOURCE_DIR}/Graph/graph.cpp)
add_executable(Sort ${PROJECT_SOURCE_DIR}/Sort/sortsAlgorithms.cpp)
