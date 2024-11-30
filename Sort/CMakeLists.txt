# 设置项目名称
cmake_minimum_required(VERSION 3.10)
project(SortProject)

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# 设置源文件
set(SOURCE_FILES
    sortsAlgorithms.cpp
)

# 添加可执行文件
add_executable(Sort ${SOURCE_FILES})
