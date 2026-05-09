#include "heat_solver.hpp"

#include <vector>
#include <string>
#include <iostream>
#include <utility>
#include <chrono>

int main(int argc, char **argv)
{
    int nx = 64;
    int ny = 64;
    double alpha = 0.25;
    int num_steps = 500;

    if (argc >= 2)
    {
        nx = std::stoi(argv[1]);
        ny = nx;
    }
    if (argc >= 3)
    {
        num_steps = std::stoi(argv[2]);
    }
    if (argc >= 4)
    {
        if (std::stod(argv[3]) > 0.25)
        {
            std::cerr << "alpha must be less than or equal to 0.25" << "\n";
        }
        alpha = std::stod(argv[3]);
    }
    std::vector<double> u_new(nx * ny, 0.0);
    std::vector<double> u_old(nx * ny, 0.0);

    // 线性部分数据
    initial_t_field(u_old, nx, ny);
    initial_t_field(u_new, nx, ny);
    // write_csv(u_old, "results/heat_initial.csv", nx, ny);
    auto start = std::chrono::high_resolution_clock::now();
    for (int step = 0; step < num_steps; step++)
    {
        step_heat_serial(u_old, u_new, nx, ny, alpha);
        std::swap(u_old, u_new);
    }
    auto end = std::chrono::high_resolution_clock::now();
    double runtime = std::chrono::duration<double>(end - start).count();
    double MLUPS = (num_steps * (nx - 2) * (ny - 2)) / (runtime * 1000000);

    write_csv(u_old, "results/heat_final_serial.csv", nx, ny);
    std::cout << "serial results:";
    std::cout << "Grid: " << nx << " x " << ny << "\n";
    std::cout << "Steps: " << num_steps << "\n";
    std::cout << "Alpha: " << alpha << "\n";
    std::cout << "Run_time: " << runtime << "\n";
    std::cout << "MLUPS: " << MLUPS << "\n";
    // std::cout << "Wrote results/heat_final_serial.csv\n";
    print_field_stats(u_old);
    std::cout << "\n";

    // 并行部分数据
    initial_t_field(u_old, nx, ny);
    initial_t_field(u_new, nx, ny);
    // write_csv(u_old, "results/heat_initial.csv", nx, ny);
    start = std::chrono::high_resolution_clock::now();
    for (int step = 0; step < num_steps; step++)
    {
        step_heat_omp(u_old, u_new, nx, ny, alpha);
        std::swap(u_old, u_new);
    }
    end = std::chrono::high_resolution_clock::now();
    runtime = std::chrono::duration<double>(end - start).count();
    MLUPS = (num_steps * (nx - 2) * (ny - 2)) / (runtime * 1000000);

    write_csv(u_old, "results/heat_final_omp.csv", nx, ny);
    std::cout << "parallel results:";
    std::cout << "Grid: " << nx << " x " << ny << "\n";
    std::cout << "Steps: " << num_steps << "\n";
    std::cout << "Alpha: " << alpha << "\n";
    std::cout << "Run_time: " << runtime << "\n";
    std::cout << "MLUPS: " << MLUPS << "\n";
    print_field_stats(u_old);

    return 0;
}
