#include "poisson_solver.hpp"
#include <vector>
#include <string>
#include <iostream>
#include <utility>
#ifdef _OPENMP
#include <omp.h>
#endif

int main()
{
    int nx = 200;
    int ny = 200;
    std::vector<double> u_old(nx * ny, 0.0);
    std::vector<double> u_new(nx * ny, 0.0);

    initial_t_field(u_old, 100, 10, 20, 40, nx, ny);
    initial_t_field(u_new, 100, 10, 20, 40, nx, ny);

    // write_field(u_old, "results/poisson_start.csv", nx, ny);
    // jacobi_update_field(u_old, u_new, nx, ny, 1e-6);
    // write_field(u_new, "results/poisson_final.csv", nx, ny);

    // gauss_update_field(u_old, nx, ny, 1e-6);
    SolveResult gauss_results_serial = solve_gauss(u_old, nx, ny, 1e-7);
    std::cout << "serial run time :" << gauss_results_serial.runtime << "\n";

    std::vector<int> threads = {1, 2, 4, 8, 16};
    SolveResult gauss_results_omp;
#ifdef _OPENMP
    omp_set_dynamic(0);
#endif
    for (int thread : threads)
    {
        initial_t_field(u_old, 100, 10, 20, 40, nx, ny);
#ifdef _OPENMP
        omp_set_num_threads(thread);
        gauss_results_omp = solve_gauss_omp(u_old, nx, ny, 1e-7);
#endif
        std::cout << thread << " threads omp run time :" << gauss_results_omp.runtime
                  << "  " << gauss_results_serial.runtime / gauss_results_omp.runtime << "\n";
    }

    // SolveResult gauss_results = solve_gauss(u_old, nx, ny, 1e-6);
    // initial_t_field(u_old, 100, 10, 20, 40, nx, ny);
    // SolveResult SOR_results = solve_SOR(u_old, nx, ny, 1e-6, 1.5);

    // std::cout << jacobi_results.converged << "," << jacobi_results.final_diff << "," << jacobi_results.iterations << ",";
    // std::cout << jacobi_results.mlups << "," << jacobi_results.runtime << "\n";
    // std::cout << gauss_results.converged << "," << gauss_results.final_diff << "," << gauss_results.iterations << ",";
    // std::cout << gauss_results.mlups << "," << gauss_results.runtime << "\n";
    // std::cout << SOR_results.converged << "," << SOR_results.final_diff << "," << SOR_results.iterations << ",";
    // std::cout << SOR_results.mlups << "," << SOR_results.runtime << "\n";
}
