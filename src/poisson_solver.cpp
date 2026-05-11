#include "poisson_solver.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <utility>
#include <cmath>
#include <chrono>

int idx(int x, int y, int nx)
{
    // 输入矩阵中一个元素的位置，和矩阵的x长度，得到存储在vector中的位置
    return (y * nx + x);
}

void initial_t_field(std::vector<double> &u,
                     double t_top,
                     double t_bottom,
                     double t_left,
                     double t_right,
                     int nx,
                     int ny)
{
    // 初始化一个温度场，需要输入一个初始化为0的vector，以及四个边界温度，直接引用目标vector，四个角设置上下平面优先
    for (int x = 0; x < nx; x++)
    {
        u[idx(x, 0, nx)] = t_bottom;
        u[idx(x, ny - 1, nx)] = t_top;
    }
    for (int y = 1; y < ny - 1; y++)
    {
        u[idx(0, y, nx)] = t_left;
        u[idx(nx - 1, y, nx)] = t_right;
    }
    for (int y = 1; y < ny - 1; y++)
    {
        for (int x = 1; x < nx; x++)
        {
            u[idx(x, y, nx)] = 0;
        }
    }
}

void write_field(const std::vector<double> &u, const std::string &filename, int nx, int ny)
{
    // 输入一个温度场，文件名，以及网格长度宽度，将温度场写入到文件名中
    std::ofstream file(filename);
    for (int y = 0; y < ny; y++)
    {
        for (int x = 0; x < nx; x++)
        {
            file << u[idx(x, y, nx)];
            if (x != nx - 1)
            {
                file << ",";
            }
            else
            {
                file << "\n";
            }
        }
    }
}

SolveResult solve_jacobi(std::vector<double> &u_old,
                         std::vector<double> &u_new,
                         int nx,
                         int ny,
                         double dif_stop)
{
    // 用来按照中心点温度是四周的平均值方式更新温度,直至最大误差小于目标误差,并且记录一个迭代误差的过程。
    double dif_now = 0;
    int i = 0;
    // std::ofstream file("results/jacobi_convence.csv");
    auto start = std::chrono::high_resolution_clock::now();
    do
    {
        dif_now = 0;
        i = i + 1;
        for (int y = 1; y < ny - 1; y++)
        {
            for (int x = 1; x < nx - 1; x++)
            {
                u_new[idx(x, y, nx)] = (u_old[idx(x - 1, y, nx)] +
                                        u_old[idx(x + 1, y, nx)] +
                                        u_old[idx(x, y - 1, nx)] +
                                        u_old[idx(x, y + 1, nx)]) /
                                       4;
                if (std::abs(u_new[idx(x, y, nx)] - u_old[idx(x, y, nx)]) > dif_now)
                {
                    dif_now = std::abs(u_new[idx(x, y, nx)] - u_old[idx(x, y, nx)]);
                }
            }
        }
        std::swap(u_old, u_new);
        // file << i << "," << dif_now << "\n";
        // std::cout << i << "," << dif_now << "\n";

    } while ((dif_now > dif_stop) && i < 50000);
    auto end = std::chrono::high_resolution_clock::now();
    double runtime = std::chrono::duration<double>(end - start).count();
    bool converged = true;
    double mlups = (i * (nx - 2) * (ny - 2)) / (1000000 * runtime);
    if (i == 50000)
    {
        converged = false;
    }
    return SolveResult{
        i,
        dif_now,
        runtime,
        mlups,
        converged};
}

SolveResult solve_gauss(std::vector<double> &u,
                        int nx,
                        int ny,
                        double dif_stop)
{
    // 用来按照中心点温度是四周的平均值方式更新温度,直至最大误差小于目标误差,并且记录一个迭代误差的过程。
    double dif_now = 0;
    int i = 0;
    std::ofstream file("results/gauss_convence.csv");
    auto start = std::chrono::high_resolution_clock::now();
    do
    {
        dif_now = 0;
        i = i + 1;
        for (int y = 1; y < ny - 1; y++)
        {
            for (int x = 1; x < nx - 1; x++)
            {
                double old_value = u[idx(x, y, nx)];
                u[idx(x, y, nx)] = (u[idx(x - 1, y, nx)] +
                                    u[idx(x + 1, y, nx)] +
                                    u[idx(x, y - 1, nx)] +
                                    u[idx(x, y + 1, nx)]) /
                                   4;
                if (std::abs(old_value - u[idx(x, y, nx)]) > dif_now)
                {
                    dif_now = std::abs(old_value - u[idx(x, y, nx)]);
                }
            }
        }
        file << i << "," << dif_now << "\n";
        // std::cout << i << "," << dif_now << "\n";

    } while ((dif_now > dif_stop) && i < 50000);
    auto end = std::chrono::high_resolution_clock::now();
    double runtime = std::chrono::duration<double>(end - start).count();
    bool converged = true;
    double mlups = (i * (nx - 2) * (ny - 2)) / (1000000 * runtime);
    if (i == 50000)
    {
        converged = false;
    }
    return SolveResult{
        i,
        dif_now,
        runtime,
        mlups,
        converged};
}

SolveResult solve_SOR(std::vector<double> &u,
                      int nx,
                      int ny,
                      double dif_stop,
                      double omega)
{
    // 用来按照中心点温度是四周的平均值方式更新温度,直至最大误差小于目标误差,并且记录一个迭代误差的过程。
    double dif_now = 0;
    int i = 0;
    std::ofstream file("results/SOR_convence.csv");
    auto start = std::chrono::high_resolution_clock::now();
    do
    {
        dif_now = 0;
        i = i + 1;
        for (int y = 1; y < ny - 1; y++)
        {
            for (int x = 1; x < nx - 1; x++)
            {
                double old_value = u[idx(x, y, nx)];
                double dif = (u[idx(x - 1, y, nx)] +
                              u[idx(x + 1, y, nx)] +
                              u[idx(x, y - 1, nx)] +
                              u[idx(x, y + 1, nx)]) /
                                 4 -
                             old_value;
                u[idx(x, y, nx)] = omega * dif + old_value;
                if (std::abs(old_value - u[idx(x, y, nx)]) > dif_now)
                {
                    dif_now = std::abs(old_value - u[idx(x, y, nx)]);
                }
            }
        }
        file << i << "," << dif_now << "\n";
        // std::cout << i << "," << dif_now << "\n";

    } while ((dif_now > dif_stop) && i < 50000);
    auto end = std::chrono::high_resolution_clock::now();
    double runtime = std::chrono::duration<double>(end - start).count();
    bool converged = true;
    double mlups = (i * (nx - 2) * (ny - 2)) / (1000000 * runtime);
    if (i == 50000)
    {
        converged = false;
    }
    return SolveResult{
        i,
        dif_now,
        runtime,
        mlups,
        converged};
}

SolveResult solve_jacobi_omp(std::vector<double> &u_old,
                             std::vector<double> &u_new,
                             int nx,
                             int ny,
                             double dif_stop)
{
    // 用来按照中心点温度是四周的平均值方式更新温度,直至最大误差小于目标误差,并且记录一个迭代误差的过程。
    double dif_now = 0;
    int i = 0;
    std::ofstream file("results/jacobi_convence.csv");
    auto start = std::chrono::high_resolution_clock::now();
    do
    {
        dif_now = 0;
        i = i + 1;
#pragma omp parallel for reduction(max : dif_now) schedule(static)
        for (int y = 1; y < ny - 1; y++)
        {
            for (int x = 1; x < nx - 1; x++)
            {
                u_new[idx(x, y, nx)] = (u_old[idx(x - 1, y, nx)] +
                                        u_old[idx(x + 1, y, nx)] +
                                        u_old[idx(x, y - 1, nx)] +
                                        u_old[idx(x, y + 1, nx)]) /
                                       4;
                if (std::abs(u_new[idx(x, y, nx)] - u_old[idx(x, y, nx)]) > dif_now)
                {
                    dif_now = std::abs(u_new[idx(x, y, nx)] - u_old[idx(x, y, nx)]);
                }
            }
        }
        std::swap(u_old, u_new);
        // file << i << "," << dif_now << "\n";
        // std::cout << i << "," << dif_now << "\n";

    } while ((dif_now > dif_stop) && i < 50000);
    auto end = std::chrono::high_resolution_clock::now();
    double runtime = std::chrono::duration<double>(end - start).count();
    bool converged = true;
    double mlups = (i * (nx - 2) * (ny - 2)) / (1000000 * runtime);
    if (i == 50000)
    {
        converged = false;
    }
    return SolveResult{
        i,
        dif_now,
        runtime,
        mlups,
        converged};
}

SolveResult solve_gauss_omp(std::vector<double> &u,
                            int nx,
                            int ny,
                            double dif_stop)
{
    // 用来按照中心点温度是四周的平均值方式更新温度,直至最大误差小于目标误差,并且记录一个迭代误差的过程。
    double dif_now = 0;
    int i = 0;
    // std::ofstream file("results/gauss_convence.csv");
    auto start = std::chrono::high_resolution_clock::now();
    do
    {
        dif_now = 0;
        i = i + 1;
#pragma omp parallel for reduction(max : dif_now) schedule(static)
        for (int y = 1; y < ny - 1; y++)
        {
            for (int x = 1; x < nx - 1; x++)
            {
                if ((x + y) % 2 == 0)
                {
                    double old_value = u[idx(x, y, nx)];
                    u[idx(x, y, nx)] = (u[idx(x - 1, y, nx)] +
                                        u[idx(x + 1, y, nx)] +
                                        u[idx(x, y - 1, nx)] +
                                        u[idx(x, y + 1, nx)]) /
                                       4;
                    if (std::abs(old_value - u[idx(x, y, nx)]) > dif_now)
                    {
                        dif_now = std::abs(old_value - u[idx(x, y, nx)]);
                    }
                }
            }
        }
#pragma omp parallel for reduction(max : dif_now) schedule(static)
        for (int y = 1; y < ny - 1; y++)
        {
            for (int x = 1; x < nx - 1; x++)
            {
                if ((x + y) % 2 == 1)
                {
                    double old_value = u[idx(x, y, nx)];
                    u[idx(x, y, nx)] = (u[idx(x - 1, y, nx)] +
                                        u[idx(x + 1, y, nx)] +
                                        u[idx(x, y - 1, nx)] +
                                        u[idx(x, y + 1, nx)]) /
                                       4;
                    if (std::abs(old_value - u[idx(x, y, nx)]) > dif_now)
                    {
                        dif_now = std::abs(old_value - u[idx(x, y, nx)]);
                    }
                }
            }
        }

        // file << i << "," << dif_now << "\n";
        // std::cout << i << "," << dif_now << "\n";

    } while ((dif_now > dif_stop) && i < 50000);

    auto end = std::chrono::high_resolution_clock::now();
    double runtime = std::chrono::duration<double>(end - start).count();
    bool converged = true;
    double mlups = (i * (nx - 2) * (ny - 2)) / (1000000 * runtime);
    if (i == 50000)
    {
        converged = false;
    }
    return SolveResult{
        i,
        dif_now,
        runtime,
        mlups,
        converged};
}
