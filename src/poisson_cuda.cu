// 1. include
#include "poisson_solver.hpp"
#include <vector>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <chrono>
#include <cuda_runtime.h>
#include <string>

{
    std::cout << __PRETTY_FUNCTION__ << "\n";
}

// 2. cuda kernels
__global__ void update_jacobi(float *u_old,
                              float *u_new,
                              float *diff,
                              int nx,
                              int ny)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    if (x > 0 && x < nx - 1 && y > 0 && y < ny - 1)
    {
        u_new[idx(x, y, nx)] = (u_old[idx(x - 1, y, nx)] +
                                u_old[idx(x + 1, y, nx)] +
                                u_old[idx(x, y - 1, nx)] +
                                u_old[idx(x, y + 1, nx)]) /
                               4;
        diff[idx(x, y, nx)] = fabsf(u_new[idx(x, y, nx)] - u_old[idx(x, y, nx)]);
    }
}

__global__ void update_diff_block_2d(const float *diff,
                                     float *diff_block,
                                     int nx,
                                     int ny)
{
    extern __shared__ float sdata[];
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int local_id = threadIdx.y * blockDim.x + threadIdx.x;
    int threads_per_block = blockDim.x * blockDim.y;
    float value = 0.0f;
    if (x < nx && y < ny)
    {
        value = diff[idx(x, y, nx)];
    }
    sdata[local_id] = value;
    __syncthreads();
    for (int stride = threads_per_block / 2; stride > 0; stride /= 2)
    {
        if (local_id < stride)
        {
            sdata[local_id] = fmaxf(sdata[local_id], sdata[local_id + stride]);
        }
        __syncthreads();
    }
    if (local_id == 0)
    {
        int block_id = blockIdx.y * gridDim.x + blockIdx.x;
        diff_block[block_id] = sdata[0];
    }
}

__global__ void reduce_diff(const float *diff1,
                            float *diff2,
                            int dif_dim)
{
    extern __shared__ float sdata[];
    int local_id = threadIdx.x;
    int num = blockIdx.x * blockDim.x + threadIdx.x;
    float value = 0.0f;
    if (num < dif_dim)
    {
        value = diff1[num];
    }
    sdata[local_id] = value;
    __syncthreads();
    for (int stride = blockDim.x / 2; stride > 0; stride /= 2)
    {
        if (local_id < stride)
        {
            sdata[local_id] = fmaxf(sdata[local_id], sdata[local_id + stride]);
        }
        __syncthreads();
    }
    if (local_id == 0)
    {
        diff2[blockIdx.x] = sdata[0];
    }
}

// 3. host helper function
float reduce_max_dif(float *buffer_diff1,
                     float *buffer_diff2,
                     int n,
                     int block_size)
{
    int current_n = n;
    while (current_n > 1)
    {
        int grid_size = (current_n + block_size - 1) / block_size;
        reduce_diff<<<grid_size, block_size, block_size * sizeof(float)>>>(buffer_diff1, buffer_diff2, current_n);

        current_n = grid_size;
        std::swap(buffer_diff1, buffer_diff2);
    }
    float max_diff = 0.0f;
    cudaMemcpy(&max_diff, buffer_diff1, sizeof(float), cudaMemcpyDeviceToHost);
    return max_diff;
}

void jacobi_cuda(float *&d_u_old,
                 float *&d_u_new,
                 float *d_diff,
                 std::vector<double> &u_old,
                 std::vector<float> &u_old_f,
                 dim3 block_dim,
                 dim3 grid_dim,
                 float goal_diff,
                 int block_size,
                 int write,
                 float *d_buffer_diff1,
                 float *d_buffer_diff2,
                 int &iteration,
                 float &max_diff_now,
                 int max_iteration,
                 int nx,
                 int ny,
                 std::size_t bytes)
{
    do
    {
        update_jacobi<<<grid_dim, block_dim>>>(d_u_old, d_u_new, d_diff, nx, ny);
        iteration += 1;
        std::swap(d_u_old, d_u_new);
        cudaMemcpy(d_buffer_diff1, d_diff, sizeof(float) * static_cast<std::size_t>(nx * ny), cudaMemcpyDeviceToDevice);

        max_diff_now = reduce_max_dif(d_buffer_diff1, d_buffer_diff2, nx * ny, block_size);
        std::cout << max_diff_now << ",";
    } while (max_diff_now > goal_diff && iteration < max_iteration);

    cudaMemcpy(u_old_f.data(), d_u_old, bytes, cudaMemcpyDeviceToHost);
    vector_float2double(u_old_f, u_old, nx, ny);

    if (write == 1)
    {
        write_field(u_old, "./results/poisson/jacobi.csv", nx, ny);
    }
}
// 4. main

int main(int argc, char **argv)
{
    // 开始计时
    auto start = std::chrono::high_resolution_clock::now();

    // 全局参数部分
    int nx = 100;
    int ny = 100;
    double temp_top = 100;
    double temp_bottom = 10;
    double temp_left = 20;
    double temp_right = 40;
    float goal_diff = 1e-4f;
    int write = 0;
    int max_iteration = 500000;

    if (argc > 1)
    {
        nx = std::stoi(argv[1]);
        ny = nx;
    }
    if (argc > 2)
    {
        write = std::stoi(argv[2]);
    }
    if (argc > 3)
    {
        goal_diff = std::stof(argv[3]);
    }
    // 初始化容器
    std::vector<double> u_old(nx * ny, 0.0);
    std::vector<double> u_new(nx * ny, 0.0);
    initial_t_field(u_old, temp_top, temp_bottom, temp_left, temp_right, nx, ny);
    initial_t_field(u_new, temp_top, temp_bottom, temp_left, temp_right, nx, ny);

    std::vector<float> u_old_f(nx * ny, 0.0f);
    std::vector<float> u_new_f(nx * ny, 0.0f);
    vector_double2float(u_old, u_old_f, nx, ny);
    vector_double2float(u_new, u_new_f, nx, ny);
    std::vector<float> h_diff(nx * ny, 0.0f);

    std::size_t bytes = static_cast<std::size_t>(nx * ny) * sizeof(float);

    float *d_u_old = nullptr;
    float *d_u_new = nullptr;
    float *d_diff = nullptr;

    cudaMalloc(&d_u_old, bytes);
    cudaMalloc(&d_u_new, bytes);
    cudaMalloc(&d_diff, bytes);

    cudaMemcpy(d_u_old, u_old_f.data(), bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_u_new, u_new_f.data(), bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_diff, h_diff.data(), bytes, cudaMemcpyHostToDevice);

    dim3 block_dim(16, 16);
    dim3 grid_dim(
        (nx + block_dim.x - 1) / block_dim.x,
        (ny + block_dim.y - 1) / block_dim.y);

    int block_size = 256;

    // 创造两个buffer diff
    float *d_buffer_diff1 = nullptr;
    float *d_buffer_diff2 = nullptr;
    cudaMalloc(&d_buffer_diff1, static_cast<std::size_t>(nx * ny) * sizeof(float));
    cudaMalloc(&d_buffer_diff2, static_cast<std::size_t>(nx * ny) * sizeof(float));

    int iteration = 0;
    float max_diff_now;

    // 执行jacobi迭代法
    jacobi_cuda(d_u_old, d_u_new, d_diff, u_old, u_old_f, block_dim, grid_dim, goal_diff, block_size, write, d_buffer_diff1, d_buffer_diff2, iteration, max_diff_now, max_iteration, nx, ny, bytes);

    auto end = std::chrono::high_resolution_clock::now();
    auto total_time_ms = std::chrono::duration<double, std::milli>(end - start).count();
    std::cout << "grid size: " << nx << " x " << ny << "\n";
    std::cout << "goal diff: " << goal_diff << "\n";
    std::cout << "final diff: " << max_diff_now << "\n";
    std::cout << "iteration: " << iteration << "\n";
    std::cout << "total time:" << total_time_ms << "ms\n";
}
