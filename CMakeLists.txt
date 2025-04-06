cmake_minimum_required(VERSION 3.10)
project(gauss_newton)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(PROJECT_SOURCE_DIR "${CMAKE_SOURCE_DIR}")
add_compile_definitions(PROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}")

add_executable(${PROJECT_NAME} src/main.cpp)
