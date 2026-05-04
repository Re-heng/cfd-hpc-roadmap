#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <utility>
#include <chrono>
#include <algorithm>
#include <numeric>

int idx(int x, int y, int nx)
{
    return nx * y + x;
}

void write_csv(const std::vector<double> &u,
               const std::string &filename,
               int nx,
               int ny)
{
    std::ofstream file(filename);
    for (int j = 0; j < ny; j++)
    {
        for (int i = 0; i < nx; i++)
        {
            if (i == nx - 1)
            {
                file << u[idx(i, j, nx)];
            }
            else
            {
                file << u[idx(i, j, nx)];
                file << ',';
            }
        }
        if (j != ny - 1)
        {
            file << '\n';
        }
    }
}

void initial_t_field(std::vector<double> &u,
                     int nx,
                     int ny)
{
    for (int j = 0; j < ny; j++)
    {
        for (int i = 0; i < nx; i++)
        {

            if (i >= 1.0 / 4.0 * nx && i < 3.0 / 4.0 * nx && j >= 1.0 / 4.0 * ny && j < 3.0 / 4.0 * ny)
            {
                u[idx(i, j, nx)] = 40;
            }
        }
    }
}

void print_field(const std::vector<double> &u,
                 int nx,
                 int ny)
{
    for (int j = 0; j <= ny - 1; j++)
    {
        for (int i = 0; i <= nx - 1; i++)
        {
            std::cout << u[idx(i, j, nx)] << ' ';
        }
        std::cout << std::endl;
    }
}
void step_heat(const std::vector<double> &u_old,
               std::vector<double> &u_new,
               int nx,
               int ny,
               double alpha)
{
    for (int j = 1; j < ny - 1; j++)
    {
        for (int i = 1; i < nx - 1; i++)
        {
            u_new[idx(i, j, nx)] = u_old[idx(i, j, nx)] + alpha * (u_old[idx(i - 1, j, nx)] + u_old[idx(i + 1, j, nx)] + u_old[idx(i, j - 1, nx)] + u_old[idx(i, j + 1, nx)] - 4 * u_old[idx(i, j, nx)]);
        }
    }
}

void print_field_stats(const std::vector<double> &u)
{
    auto min_max = std::minmax_element(u.begin(), u.end());
    double sum = std::accumulate(u.begin(), u.end(), 0.0);
    double mean = sum / static_cast<double>(u.size());
    std::cout << "min temperature: " << *min_max.first << "\n";
    std::cout << "max temperature: " << *min_max.second << "\n";
    std::cout << "mean temperature: " << mean << "\n";
}

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
    initial_t_field(u_old, nx, ny);
    write_csv(u_old, "results/heat_initial.csv", nx, ny);
    auto start = std::chrono::high_resolution_clock::now();
    for (int step = 0; step < num_steps; step++)
    {
        step_heat(u_old, u_new, nx, ny, alpha);
        std::swap(u_old, u_new);
    }
    auto end = std::chrono::high_resolution_clock::now();
    double runtime = std::chrono::duration<double>(end - start).count();
    double MLUPS = (num_steps * (nx - 2) * (ny - 2)) / (runtime * 1000000);

    write_csv(u_old, "results/heat_final.csv", nx, ny);
    std::cout << "Grid: " << nx << " x " << ny << "\n";
    std::cout << "Steps: " << num_steps << "\n";
    std::cout << "Alpha: " << alpha << "\n";
    std::cout << "Run_time: " << runtime << "\n";
    std::cout << "MLUPS: " << MLUPS << "\n";
    std::cout << "Wrote results/heat_final.csv\n";
    print_field_stats(u_old);
    return 0;
}
