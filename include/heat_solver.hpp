#pragma once

#include <string>
#include <vector>

#ifdef __CUDACC__
#define HD __host__ __device__
#define FINLINE __forceinline__
#else
#define HD
#define FINLINE inline
#endif

HD FINLINE int idx(int x, int y, int nx)
{
    return nx * y + x;
}

void write_csv(const std::vector<double> &u,
               const std::string &filename,
               int nx,
               int ny);

void initial_t_field(std::vector<double> &u,
                     int nx,
                     int ny);

void print_field(const std::vector<double> &u,
                 int nx,
                 int ny);

void step_heat_serial(const std::vector<double> &u_old,
                      std::vector<double> &u_new,
                      int nx,
                      int ny,
                      double alpha);

void step_heat_omp(const std::vector<double> &u_old,
                   std::vector<double> &u_new,
                   int nx,
                   int ny,
                   double alpha);

void print_field_stats(const std::vector<double> &u);
