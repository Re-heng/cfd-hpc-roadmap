#include <iostream>
#include <vector>

int index(int x, int y, int nx)
{
    return nx * y + x;
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
    return 0;
}
