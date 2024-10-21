#include <stdio.h>
#include <cuda_runtime.h>
#include <cmath> // For math functions
#include "pt.cuh"


// Example kernel that computes the CDF for an array of t-values
__global__ void compute_cdf(double *t_values, double *df_values, double *results, int n) {
    int idx = threadIdx.x + blockDim.x * blockIdx.x;
    if (idx < n) {
        double t = t_values[idx];
        double df = df_values[idx];

        // Call the pt function
        double cdf = pt(t, df);

        // Store the result
        results[idx] = cdf;
    }
}

int main() {
    // Number of elements
    int n = 10;

    // Host arrays
    double *h_t_values = (double*)malloc(n * sizeof(double));
    double *h_df_values = (double*)malloc(n * sizeof(double));
    double *h_results = (double*)malloc(n * sizeof(double));

    // Initialize t-values and df-values
    for (int i = 0; i < n; i++) {
        h_t_values[i] = -3.0 + i * 0.7; // Example t-values from -3 to 3
        h_df_values[i] = 10.0;          // Degrees of freedom
    }

    // Device arrays
    double *d_t_values, *d_df_values, *d_results;
    cudaMalloc((void**)&d_t_values, n * sizeof(double));
    cudaMalloc((void**)&d_df_values, n * sizeof(double));
    cudaMalloc((void**)&d_results, n * sizeof(double));

    // Copy data to device
    cudaMemcpy(d_t_values, h_t_values, n * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_df_values, h_df_values, n * sizeof(double), cudaMemcpyHostToDevice);

    // Launch kernel
    int blockSize = 256;
    int gridSize = (n + blockSize - 1) / blockSize;
    compute_cdf<<<gridSize, blockSize>>>(d_t_values, d_df_values, d_results, n);

    // Copy results back to host
    cudaMemcpy(h_results, d_results, n * sizeof(double), cudaMemcpyDeviceToHost);

    // Output the results
    for (int i = 0; i < n; i++) {
        printf("t: %f, df: %f, CDF: %f\n", h_t_values[i], h_df_values[i], (1-h_results[i]));
    }

    // Free memory
    cudaFree(d_t_values);
    cudaFree(d_df_values);
    cudaFree(d_results);
    free(h_t_values);
    free(h_df_values);
    free(h_results);

    return 0;
}
