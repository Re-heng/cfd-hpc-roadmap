#include "poisson_solver.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <utility>
#include <cmath>

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

void update_field_single(std::vector<double> &u_old,
                         std::vector<double> &u_new,
                         int nx,
                         int ny)
{
    // 用来按照中心点温度是四周的平均值方式更新温度,仅仅更新一次。
    for (int y = 1; y < ny - 1; y++)
    {
        for (int x = 1; x < nx - 1; x++)
        {
            u_new[idx(x, y, nx)] = (u_old[idx(x - 1, y, nx)] +
                                    u_old[idx(x + 1, y, nx)] +
                                    u_old[idx(x, y - 1, nx)] +
                                    u_old[idx(x, y + 1, nx)]) /
                                   4;
        }
    }
}

void update_field_final(std::vector<double> &u_old,
                        std::vector<double> &u_new,
                        int nx,
                        int ny,
                        double stop_dif)
{
    // 用来按照中心点温度是四周的平均值方式更新温度,直至更新前后对应点的最大温差绝对值小于 stop_dif
    update_field_single(u_old, u_new, nx, ny);
    double max_dif = 0;
    for (int y = 1; y < ny - 1; y++)
    {
        for (int x = 1; x < nx - 1; x++)
        {
            if (std::abs(u_new[idx(x, y, nx)] - u_old[idx(x, y, nx)]) > max_dif)
            {
                max_dif = std::abs(u_new[idx(x, y, nx)] - u_old[idx(x, y, nx)]);
            }
        }
    }
    while (max_dif > stop_dif)
    {
        std::swap(u_old, u_new);
        std::cout << "dif now :" << max_dif;
        update_field_single(u_old, u_new, nx, ny);
        max_dif = 0;
        for (int y = 1; y < ny - 1; y++)
        {
            for (int x = 1; x < nx - 1; x++)
            {
                if (std::abs(u_new[idx(x, y, nx)] - u_old[idx(x, y, nx)]) > max_dif)
                {
                    max_dif = std::abs(u_new[idx(x, y, nx)] - u_old[idx(x, y, nx)]);
                }
            }
        }
    }
}
