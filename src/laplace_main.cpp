#include "poisson_solver.hpp"
#include <vector>
#include <string>
#include <iostream>
#include <utility>

int main()
{
    int nx = 100;
    int ny = 100;
    std::vector<double> u_old(nx * ny, 0.0);
    std::vector<double> u_new(nx * ny, 0.0);

    initial_t_field(u_old, 100, 10, 20, 40, nx, ny);
    initial_t_field(u_new, 100, 10, 20, 40, nx, ny);

    // write_field(u_old, "results/poisson_start.csv", nx, ny);
    // jacobi_update_field(u_old, u_new, nx, ny, 1e-6);
    // write_field(u_new, "results/poisson_final.csv", nx, ny);

    // gauss_update_field(u_old, nx, ny, 1e-6);
    SolveResult jacobi_results = solve_jacobi(u_old, u_new, nx, ny, 1e-6);
    SolveResult gauss_results = solve_gauss(u_old, nx, ny, 1e-6);
    SolveResult SOR_results = solve_SOR(u_old, nx, ny, 1e-6, 1.5);
}
