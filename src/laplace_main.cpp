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
    write_field(u_old, "results/poisson_start.csv", nx, ny);
    update_field(u_old, u_new, nx, ny, 1e-6);
    write_field(u_new, "results/poisson_final.csv", nx, ny);
    /*   for (int i = 0; i < 10; i++)
       {·
           std::swap(u_old, u_new);
           update_field_single(u_old, u_new, nx, ny);
       }
       write_field(u_new, "results/poisson_step10.csv", nx, ny);
   */
}
