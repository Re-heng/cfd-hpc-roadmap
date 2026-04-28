#include <iostream>
#include <vector>
#include <fstream>
#include <string>

int index(int x, int y, int nx)
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
                file << u[index(i, j, nx)];
            }
            else
            {
                file << u[index(i, j, nx)];
                file << ',';
            }
        }
        if (j != ny - 1)
        {
            file << '\n';
        }
    }
}

int main()
{
    int nx = 5;
    int ny = 4;
    std::vector<double> u(nx * ny, 0.0);

    for (int j = 0; j <= ny - 1; j++)
    {
        for (int i = 0; i <= nx - 1; i++)
        {
            u[index(i, j, nx)] = i + j * 10;
        }
    }

    for (int j = 0; j <= ny - 1; j++)
    {
        for (int i = 0; i <= nx - 1; i++)
        {
            std::cout << u[index(i, j, nx)] << ' ';
        }
        std::cout << std::endl;
    }

    write_csv(u, "results/grid_demo.csv", nx, ny);
    std::cout << "Wrote results/grid_demo.csv\n";
    return 0;
}
