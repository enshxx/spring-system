cmake_minimum_required(VERSION 3.30)
project(gauss_newton)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(${PROJECT_NAME} src/main.cpp)
add_executable(ludecomptest src/ludcmp_test.cpp)

enable_testing()
add_test(NAME ludecomptest COMMAND ludecomptest)