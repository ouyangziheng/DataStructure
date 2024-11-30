# 设置项目名称
cmake_minimum_required(VERSION 3.10)
project(Graph)

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 14)

# 设置头文件目录
include_directories(${PROJECT_SOURCE_DIR}/include)

# 添加源文件
set(SOURCE_FILES
    graph.cpp
)

add_executable(Graph ${SOURCE_FILES})
