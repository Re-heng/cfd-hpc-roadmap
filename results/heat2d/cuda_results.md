# Heat2D CPU/OpenMP/CUDA Benchmark

## Problem setup

- nx = 1000
- ny = 1000
- steps = 1000
- alpha = 0.2
- updates = 998 * 998 * 1000 = 996,004,000

## CPU / OpenMP

| Method | Threads | Runtime (s) |   MLUPS |  Speedup |
| ------ | ------: | ----------: | ------: | -------: |
| serial |       1 |     20.2895 | 49.0896 |      1.0 |
| omp    |       1 |     20.3283 | 48.9959 | 0.998091 |
| omp    |       2 |     10.4833 | 95.0083 |  1.93541 |
| omp    |       4 |      5.7438 | 173.405 |  3.53242 |
| omp    |       8 |     3.55603 | 280.089 |  5.70567 |
| omp    |      16 |     3.48123 | 286.107 |  5.82827 |

## CUDA

| Phase          | Time (ms) |
| -------------- | --------: |
| Host to Device |   1.40384 |
| GPU solve      |   22.1458 |
| Device to Host |  0.680704 |
| Total          |   24.2303 |

## Derived CUDA performance

- CUDA kernel time = 0.0221458 s
- CUDA total time = 0.0242303 s
- CUDA kernel-only MLUPS ≈ 44,974
- CUDA total MLUPS ≈ 41,106
- CUDA total speedup vs serial ≈ 837x
- CUDA total speedup vs OpenMP 16 ≈ 144x

## Correctness validation

```text
max difference between serial and cuda is 0.00030000000000285354
max difference between serial and openmp is 0.0
max difference between cuda and openmp is 0.00030000000000285354
