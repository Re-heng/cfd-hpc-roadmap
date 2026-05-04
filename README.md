# CFD HPC Roadmap

This repository records my learning path in C++, numerical methods, CFD, and HPC.

## Project 1: 2D Heat Diffusion Solver

This project implements a simple explicit finite difference solver for the 2D heat equation.

The goal is to practice:

- C++ numerical programming
- CMake project workflow
- 2D grid indexing with a 1D `std::vector<double>`
- Explicit finite difference methods
- CSV output and Python visualization
- Basic performance measurement

## Features

- 2D temperature field stored in a 1D `std::vector<double>`
- 2D-to-1D index mapping with `idx(i, j, nx)`
- Center hot square initialization
- Explicit 5-point stencil heat diffusion update
- Double buffering with `u_old` and `u_new`
- Command line arguments for grid size, step count, and alpha
- Runtime timing with `std::chrono`
- MLUPS performance metric
- Min / max / mean temperature diagnostics
- Python visualization using matplotlib

## Numerical Method

The update formula is:

u_new(i,j) = u_old(i,j) + alpha * (
    u_old(i-1,j) + u_old(i+1,j)
  + u_old(i,j-1) + u_old(i,j+1)
  - 4 * u_old(i,j)
)
This is an explicit finite difference scheme using a 5-point stencil.

For the 2D explicit heat equation, the stability condition is approximately:

alpha <= 0.25

When alpha > 0.25, the simulation may become unstable and show a checkerboard-like numerical oscillation.


## Project Structure
.
├── CMakeLists.txt
├── README.md
├── src/
│   └── main.cpp
├── scripts/
│   └── plot_heat.py
├── docs/
│   └── progress_log.md
├── results/
│   ├── heat_initial.png
│   └── heat_final.png
└── environment.yml

## Requirements

C++:

CMake
A C++17 compiler

Python:

Python 3.10
NumPy
Matplotlib

The Python environment can be created using Conda:

conda env create -f environment.yml
conda activate cfd-hpc
Build
cmake -S . -B build
cmake --build build

## Run

Default run:

./build/heat2d

Custom run:

./build/heat2d <grid_size> <num_steps> <alpha>

Example:

./build/heat2d 128 1000 0.1

Arguments:

grid_size   Number of grid points in x and y directions
num_steps   Number of explicit time steps
alpha       Diffusion parameter

## Output

The solver writes:

results/heat_initial.csv
results/heat_final.csv

It also prints runtime and field diagnostics:

Grid: 128 x 128
Steps: 1000
Alpha: 0.1
Runtime: ...
MLUPS: ...
Min temperature: ...
Max temperature: ...
Mean temperature: ...

## Visualization
conda activate cfd-hpc
python scripts/plot_heat.py results/heat_initial.csv
python scripts/plot_heat.py results/heat_final.csv

The generated figures are:

results/heat_initial.png
results/heat_final.png

## Example Result

The initial condition is a hot square in the center of the domain. After diffusion, the heat spreads outward and the maximum temperature decreases.

## Stability Observation

Using a stable value such as:

./build/heat2d 128 1000 0.1

produces a smooth diffusion result.

Using an unstable value such as:

./build/heat2d 128 100 0.3

may produce checkerboard-like numerical oscillations.

This demonstrates that numerical simulations require stability constraints, not just correct-looking code.

## Current Status

This project is a minimal complete CPU-based 2D heat diffusion solver.

## Future Work
Split the code into .hpp and .cpp files
Add benchmark tables for different grid sizes
Add OpenMP parallelization
Implement a 2D Poisson solver
