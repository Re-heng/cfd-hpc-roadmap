cmake_minimum_required(VERSION 3.16)
project(cfd_hpc_raodmap LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(heat2d
    src/heat_main.cpp
    src/heat_solver.cpp
)

add_executable(poisson
    src/laplace_main.cpp
    src/poisson_solver.cpp)

target_include_directories(heat2d PRIVATE include)
target_include_directories(poisson PRIVATE include)
