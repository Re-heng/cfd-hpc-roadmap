#include "poisson_solver.hpp"
#include <vector>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <chrono>
#include <cuda_runtime.h>
#include <string>

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

int main(int argc, char **argv)
{
    // 开始计时
    auto start = std::chrono::high_resolution_clock::now();
    int nx = 100;
    int ny = 100;
    double temp_top = 100;
    double temp_bottom = 10;
    double temp_left = 20;
    double temp_right = 40;
    float goal_diff = 1e-4f;
    int write = 0;

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

    std::size_t shared_bytes = block_size * sizeof(float);

    // 创造两个buffer diff
    float *d_buffer_diff1 = nullptr;
    float *d_buffer_diff2 = nullptr;
    cudaMalloc(&d_buffer_diff1, static_cast<std::size_t>(nx * ny) * sizeof(float));
    cudaMalloc(&d_buffer_diff2, static_cast<std::size_t>(nx * ny) * sizeof(float));
    int diff_dim;

    float max_diff = 0.0f;
    int iteration = 0;

    do
    {
        update_jacobi<<<grid_dim, block_dim>>>(d_u_old, d_u_new, d_diff, nx, ny);
        iteration += 1;
        std::swap(d_u_old, d_u_new);
        cudaMemcpy(d_buffer_diff1, d_diff, sizeof(float) * static_cast<std::size_t>(nx * ny), cudaMemcpyDeviceToDevice);

        diff_dim = nx * ny;
        int reduce_grid_size = ((nx * ny) + block_size - 1) / block_size;
        while (diff_dim > 1)
        {
            reduce_diff<<<reduce_grid_size, block_size, shared_bytes>>>(d_buffer_diff1, d_buffer_diff2, diff_dim);
            diff_dim = reduce_grid_size;
            reduce_grid_size = (diff_dim + block_size - 1) / block_size;
            std::swap(d_buffer_diff1, d_buffer_diff2);
        }

        cudaMemcpy(&max_diff,
                   d_buffer_diff1,
                   sizeof(float),
                   cudaMemcpyDeviceToHost);
    } while (max_diff > goal_diff);

    cudaMemcpy(u_old_f.data(), d_u_old, bytes, cudaMemcpyDeviceToHost);
    vector_float2double(u_old_f, u_old, nx, ny);

    if (write == 1)
    {
        write_field(u_old, "./results/poisson/jacobi.csv", nx, ny);
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto total_time_ms = std::chrono::duration<double, std::milli>(end - start).count();
    std::cout << "grid size: " << nx << " x " << ny << "\n";
    std::cout << "goal diff: " << goal_diff << "\n";
    std::cout << "final diff: " << max_diff << "\n";
    std::cout << "iteration: " << iteration << "\n";
    std::cout << "total time:" << total_time_ms << "ms\n";
}
