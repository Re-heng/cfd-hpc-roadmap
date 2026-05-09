#pragma once

#include <string>
#include <vector>

int idx(int x, int y, int nx);

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
