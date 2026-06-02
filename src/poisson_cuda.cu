#include "poisson_solver.hpp"
#include <vector>
#include <fstream>
#include <algorithm>
#include <iostream>

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

int main()
{
    int nx = 500;
    int ny = 500;
    double temp_top = 100;
    double temp_bottom = 50;
    double temp_left = 40;
    double temp_right = 0;

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

    std::size_t shared_bytes =
        block_dim.x * block_dim.y * sizeof(float);

    std::vector<float> h_diff_block(grid_dim.x * grid_dim.y, 0.0f);
    float *d_diff_block = nullptr;
    cudaMalloc(&d_diff_block, static_cast<std::size_t>(grid_dim.x * grid_dim.y) * sizeof(float));

    float max_diff;
    do
    {
        update_jacobi<<<grid_dim, block_dim>>>(d_u_old, d_u_new, d_diff, nx, ny);
        std::swap(d_u_old, d_u_new);
        update_diff_block_2d<<<grid_dim, block_dim, shared_bytes>>>(d_diff, d_diff_block, nx, ny);
        cudaMemcpy(h_diff_block.data(),
                   d_diff_block,
                   static_cast<std::size_t>(grid_dim.x * grid_dim.y) * sizeof(float),
                   cudaMemcpyDeviceToHost);
        max_diff = *std::max_element(h_diff_block.begin(), h_diff_block.end());
        // std::cout << max_diff << ",";
    } while (max_diff > 1e-4);

    cudaMemcpy(u_old_f.data(), d_u_old, bytes, cudaMemcpyDeviceToHost);
    vector_float2double(u_old_f, u_old, nx, ny);
    write_field(u_old, "./results/poisson/jacobi.csv", nx, ny);
}
