## Observations

For a 1000 x 1000 grid with 1000 time steps, the OpenMP version shows clear speedup up to 8 threads.

The best result is obtained with 8 threads:
- runtime: 3.36564 s
- MLUPS: 295.933
- speedup: 6.02x

Using 16 threads is slower than using 8 threads. This is expected on the test machine, which has 8 physical cores and 16 logical threads. The 5-point heat stencil is memory-bandwidth bound, so SMT threads do not necessarily improve performance once the physical cores already saturate memory bandwidth.

The OpenMP 1-thread run is slightly slower than the serial version due to OpenMP runtime overhead.

## data
method, thread, runtime, mlups, speedup
serial,1,20.2447,49.1983,1
omp,1,21.2705,46.8255,0.951771
omp,2,10.9792,90.7172,1.84391
omp,4,5.6362,176.716,3.5919
omp,8,3.36564,295.933,6.0151
omp,16,4.10911,242.389,4.92678
