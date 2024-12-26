cmake_minimum_required(VERSION 3.10)
project(gauss_newton)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(${PROJECT_NAME} src/main.cpp)
add_executable(ludecomptest src/ludcmp_test.cpp)
add_executable(jacplusodetest src/jacplusodetest.cpp)
add_executable(qrdecomptest src/qrdecompsolver_test.cpp)
add_executable(choleskysoltest src/choleskysolver_test.cpp)
