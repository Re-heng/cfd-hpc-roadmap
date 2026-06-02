#include <cuda_runtime.h>

#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include "heat_solver.hpp"

__global__ void heat_solver_kernel_dim1(float *u_old,
                                        float *u_new,
                                        int nx,
                                        int ny,
                                        float alpha)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < nx * ny)
    {
        int y = i / nx;
        int x = i - y * nx;
        if (x != 0 && x != nx - 1 && y != 0 && y != ny - 1)
        {
            u_new[idx(x, y, nx)] = u_old[idx(x, y, nx)] + alpha * (u_old[idx(x - 1, y, nx)] +
                                                                   u_old[idx(x + 1, y, nx)] +
                                                                   u_old[idx(x, y - 1, nx)] +
                                                                   u_old[idx(x, y + 1, nx)] -
                                                                   4 * u_old[idx(x, y, nx)]);
        }
    }
}

__global__ void heat_solver_kernel_dim2(float *u_old,
                                        float *u_new,
                                        int nx,
                                        int ny,
                                        float alpha)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    if (x > 0 && x < nx - 1 && y > 0 && y < ny - 1)
    {
        u_new[idx(x, y, nx)] = u_old[idx(x, y, nx)] + alpha * (u_old[idx(x - 1, y, nx)] +
                                                               u_old[idx(x + 1, y, nx)] +
                                                               u_old[idx(x, y - 1, nx)] +
                                                               u_old[idx(x, y + 1, nx)] -
                                                               4 * u_old[idx(x, y, nx)]);
    }
}

int main(int argc, char **argv)
{
    // 加入计时
    auto cpu_setup_start = std::chrono::high_resolution_clock::now();

    // 完成一波host中的两个温度场初始化。
    int nx = 500;
    int ny = 500;
    float alpha = 0.2;
    int steps = 1000;
    int cuda_dim_mode = 2;
    int write = 0;

    if (argc > 1)
    {
        nx = std::stoi(argv[1]);
        ny = nx;
    }
    if (argc > 2)
    {
        alpha = std::stof(argv[2]);
    }
    if (argc > 3)
    {
        steps = std::stoi(argv[3]);
    }
    if (argc > 4)
    {
        write = std::stoi(argv[4]);
    }
    if (argc > 5)
    {
        cuda_dim_mode = std::stoi(argv[5]);
    }

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

    auto cpu_setup_end = std::chrono::high_resolution_clock::now();
    double cpu_setup_time = std::chrono::duration<double, std::milli>(cpu_setup_end - cpu_setup_start).count();

    cudaMalloc(&d_u_old, bytes);
    cudaMalloc(&d_u_new, bytes);

    // create timing event
    cudaEvent_t start;
    cudaEvent_t stop;

    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start);
    cudaMemcpy(d_u_old, h_u_oldf.data(), bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_u_new, h_u_newf.data(), bytes, cudaMemcpyHostToDevice);
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float h2d_ms = 0.0f;
    cudaEventElapsedTime(&h2d_ms, start, stop);

    // run gpu heat solver
    int block_size = 256;
    int grid_size = (nx * ny + block_size - 1) / block_size;

    dim3 block_dim(16, 16);
    dim3 grid_dim(
        (nx + block_dim.x - 1) / block_dim.x,
        (ny + block_dim.y - 1) / block_dim.y);

    cudaEventRecord(start);

    if (cuda_dim_mode == 1)
    {
        std::cout << "use dim1 cuda" << "\n";
        for (int step = 0; step < steps; step++)
        {
            heat_solver_kernel_dim1<<<grid_size, block_size>>>(d_u_old, d_u_new, nx, ny, alpha);
            std::swap(d_u_old, d_u_new);
        }
    }

    else
    {
        std::cout << "use dim2 cuda" << "\n";
        for (int step = 0; step < steps; step++)
        {
            heat_solver_kernel_dim2<<<grid_dim, block_dim>>>(d_u_old, d_u_new, nx, ny, alpha);
            std::swap(d_u_old, d_u_new);
        }
    }

    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float d_solve_ms = 0.0f;
    cudaEventElapsedTime(&d_solve_ms, start, stop);

    cudaEventRecord(start);
    cudaMemcpy(h_u_newf.data(), d_u_old, bytes, cudaMemcpyDeviceToHost);
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float d2h_ms = 0.0f;
    cudaEventElapsedTime(&d2h_ms, start, stop);

    for (int i = 0; i < h_u_new.size(); i++)
    {
        h_u_new[i] = static_cast<double>(h_u_newf[i]);
    }

    if (write)
    {
        write_csv(h_u_new, "results/heat2d/cuda_heat.csv", nx, ny);
    }

    auto endtime = std::chrono ::high_resolution_clock ::now();
    auto total_time_ms = std::chrono::duration<double, std::milli>(endtime - cpu_setup_start).count();

    std::cout << "grid size: " << nx << " x " << ny << "\n";
    std::cout << "steps: " << steps << "\n";
    std::cout << "alpha: " << alpha << "\n";
    std::cout << "cpu_setup: " << cpu_setup_time << "ms\n";
    std::cout << "host to device time :" << h2d_ms << "ms\n";
    std::cout << "gpu solve time :" << d_solve_ms << "ms\n";
    std::cout << "device to host time " << d2h_ms << "ms\n";
    std::cout << "cuda_total time " << h2d_ms + d_solve_ms + d2h_ms << "ms\n";
    std::cout << "total time: " << total_time_ms << "ms\n";
    std::cout << "unmeasured_overhead_ms: " << total_time_ms - cpu_setup_time - (h2d_ms + d_solve_ms + d2h_ms) << "ms\n";
}
