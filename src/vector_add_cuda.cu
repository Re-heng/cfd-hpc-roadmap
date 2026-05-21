#include <cuda_runtime.h>

#include <cmath>
#include <iostream>
#include <vector>

// This function runs on the GPU.
//
// __global__ means:
// - called from CPU code
// - executed on the GPU
//
// Each CUDA thread executes this same function.
// The thread computes its own global index i,
// then handles one element: c[i] = a[i] + b[i].
__global__ void vector_add_kernel(const float* a,
                                  const float* b,
                                  float* c,
                                  int n) {
    // blockIdx.x  : which block am I in?
    // blockDim.x  : how many threads are in each block?
    // threadIdx.x : which thread am I inside this block?
    //
    // Together they give a global 1D index.
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    // The total number of launched threads may be larger than n.
    // So we must guard against out-of-bounds access.
    if (i < n) {
        c[i] = a[i] + b[i];
    }
}

int main() {
    // Number of vector elements.
    // Start with 1 << 20 = 2^20 = 1,048,576 elements.
    int n = 1 << 20;

    // Number of bytes needed for one vector.
    std::size_t bytes = static_cast<std::size_t>(n) * sizeof(float);

    // Host memory: CPU-side arrays.
    //
    // Naming convention:
    // h_a means host a.
    // h_b means host b.
    // h_c means host c.
    std::vector<float> h_a(n);
    std::vector<float> h_b(n);
    std::vector<float> h_c(n);

    // Initialize input data on CPU.
    for (int i = 0; i < n; ++i) {
        h_a[i] = 0.5f * static_cast<float>(i);
        h_b[i] = 2.0f * static_cast<float>(i);
    }

    // Device memory: GPU-side pointers.
    //
    // Naming convention:
    // d_a means device a.
    // d_b means device b.
    // d_c means device c.
    float* d_a = nullptr;
    float* d_b = nullptr;
    float* d_c = nullptr;

    // Allocate memory on the GPU.
    cudaMalloc(&d_a, bytes);
    cudaMalloc(&d_b, bytes);
    cudaMalloc(&d_c, bytes);

    // Copy input data from CPU memory to GPU memory.
    cudaMemcpy(d_a, h_a.data(), bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, h_b.data(), bytes, cudaMemcpyHostToDevice);

    // CUDA launch configuration.
    //
    // block_size:
    //   number of threads per block
    //
    // grid_size:
    //   number of blocks
    //
    // Example:
    //   n = 1,048,576
    //   block_size = 256
    //   grid_size = 4096
    //
    // Total launched threads:
    //   grid_size * block_size
    int block_size = 256;
    int grid_size = (n + block_size - 1) / block_size;

    // Launch GPU kernel.
    //
    // This means:
    //   launch grid_size blocks
    //   each block has block_size threads
    vector_add_kernel<<<grid_size, block_size>>>(d_a, d_b, d_c, n);

    // Wait until GPU work is finished.
    //
    // Kernel launch is asynchronous:
    // CPU may continue immediately unless we synchronize.
    cudaDeviceSynchronize();

    // Copy result from GPU memory back to CPU memory.
    cudaMemcpy(h_c.data(), d_c, bytes, cudaMemcpyDeviceToHost);

    // Verify correctness on CPU.
    float max_error = 0.0f;

    for (int i = 0; i < n; ++i) {
        float expected = h_a[i] + h_b[i];
        float error = std::fabs(h_c[i] - expected);

        if (error > max_error) {
            max_error = error;
        }
    }

    // Print basic information.
    std::cout << "n = " << n << "\n";
    std::cout << "block_size = " << block_size << "\n";
    std::cout << "grid_size = " << grid_size << "\n";
    std::cout << "total_threads = " << grid_size * block_size << "\n";
    std::cout << "max_error = " << max_error << "\n";

    if (max_error < 1e-5f) {
        std::cout << "PASS\n";
    } else {
        std::cout << "FAIL\n";
    }

    // Free GPU memory.
    cudaFree(d_a);
    cudaFree(d_b);
    cudaFree(d_c);

    return 0;
}