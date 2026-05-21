cmake_minimum_required(VERSION 3.18)
project(cfd_hpc_roadmap LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# RTX 3070 Laptop GPU is Ampere, compute capability 8.6
set(CMAKE_CUDA_ARCHITECTURES 86)

find_package(OpenMP REQUIRED)

add_executable(heat2d
    src/heat_main.cpp
    src/heat_solver.cpp
)

add_executable(poisson
    src/laplace_main.cpp
    src/poisson_solver.cpp
)

add_executable(vector_add_cuda
    src/vector_add_cuda.cu
)

target_include_directories(heat2d PRIVATE include)
target_link_libraries(heat2d PRIVATE OpenMP::OpenMP_CXX)
target_include_directories(poisson PRIVATE include)
target_link_libraries(poisson PRIVATE OpenMP::OpenMP_CXX)
