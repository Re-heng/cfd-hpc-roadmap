#include <cuda_runtime.h>

#include <cmath>
#include <iostream>
#include <vector>
#include "heat_solver.hpp"

__global__ void heat_solver_kernel(float *u_old,
                                   float *u_new,
                                   int nx,
                                   int ny)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < nx * ny)
    {
        int y = i / nx;
        int x = i - y * nx;
        if (x != 0 && x != nx - 1 && y != 0 && y != ny - 1)
        {
            u_new[idx(x, y, nx)] = u_old[idx(x, y, nx)] + 0.2f * (u_old[idx(x - 1, y, nx)] +
                                                                  u_old[idx(x + 1, y, nx)] +
                                                                  u_old[idx(x, y - 1, nx)] +
                                                                  u_old[idx(x, y + 1, nx)] -
                                                                  4 * u_old[idx(x, y, nx)]);
        }
    }
}

int main()
{
    // 完成一波host中的两个温度场初始化。
    int nx = 100;
    int ny = 100;
    std::vector<double> h_u_old(nx * ny, 0.0);
    std::vector<double> h_u_new(nx * ny, 0.0);
    initial_t_field(h_u_old, nx, ny);

    std::vector<float> h_u_oldf(nx * ny, 0.0);
    std::vector<float> h_u_newf(nx * ny, 0.0);
    for (int i = 0; i < h_u_old.size(); i++)
    {
        h_u_oldf[i] = static_cast<float>(h_u_old[i]);
    }

    // 将两个温度场搬到gpu
    std::size_t bytes = static_cast<std::size_t>(nx * ny) * sizeof(float);
    float *d_u_old = nullptr;
    float *d_u_new = nullptr;

    cudaMalloc(&d_u_old, bytes);
    cudaMalloc(&d_u_new, bytes);

    cudaMemcpy(d_u_old, h_u_oldf.data(), bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_u_new, h_u_newf.data(), bytes, cudaMemcpyHostToDevice);

    // run gpu heat solver
    int block_size = 256;
    int grid_size = (nx * ny + block_size - 1) / block_size;

    for (int step = 0; step < 1000; step++)
    {
        heat_solver_kernel<<<grid_size, block_size>>>(d_u_old, d_u_new, nx, ny);
        std::swap(d_u_old, d_u_new);
    }

    cudaMemcpy(h_u_newf.data(), d_u_old, bytes, cudaMemcpyDeviceToHost);
    for (int i = 0; i < h_u_new.size(); i++)
    {
        h_u_new[i] = static_cast<double>(h_u_newf[i]);
    }
    write_csv(h_u_new, "results/cuda_heatsolver.csv", nx, ny);
}