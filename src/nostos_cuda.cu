// Nostos CUDA backend (see nostos_gpu.h). Build with --fmad=false: the CPU builds do not contract
// a*b+c into FMAs, and keeping the same rounding keeps GPU and CPU results within rounding.
#include <cstdio>
#include <cuda_runtime.h>
#include "nostos_gpu.h"

#define NG_GLOBAL
#define NG_DEVICE __device__ static inline
#define NG_KERNEL __global__
#define NG_TID ((int)(blockIdx.x * blockDim.x + threadIdx.x))
#include "nostos_kernel_body.h"

static int g_ready = 0;

#define CK(x) do { cudaError_t e_ = (x); if (e_ != cudaSuccess) { \
    fprintf(stderr, "nostos_cuda: %s failed: %s\n", #x, cudaGetErrorString(e_)); return 1; } } while (0)

extern "C" int ng_cuda_init(int device) {
    int n = 0;
    CK(cudaGetDeviceCount(&n));
    if (n <= 0) { fprintf(stderr, "nostos_cuda: no CUDA device\n"); return 1; }
    if (device < 0 || device >= n) device = 0;
    CK(cudaSetDevice(device));
    cudaDeviceProp prop;
    CK(cudaGetDeviceProperties(&prop, device));
    fprintf(stderr, "nostos_cuda: device %d: %s (sm_%d%d, %d SMs)\n", device, prop.name, prop.major, prop.minor, prop.multiProcessorCount);
    g_ready = 1;
    return 0;
}

extern "C" int ng_cuda_run(const NgBatch *b, int steps_per_launch, double *out_min_d, double *out_t_min,
                           NgProgress progress, void *user) {
    if (!g_ready) return 1;
    const int n = b->n_pairs;
    double *d_star = 0, *d_iso = 0, *d_state = 0, *d_best = 0, *d_params = 0;
    CK(cudaMalloc(&d_star, sizeof(double) * 6 * n));
    CK(cudaMalloc(&d_iso, sizeof(double) * 6 * n));
    CK(cudaMalloc(&d_state, sizeof(double) * 12 * n));
    CK(cudaMalloc(&d_best, sizeof(double) * 2 * n));
    CK(cudaMalloc(&d_params, sizeof(double) * NG_NPARAMS));
    CK(cudaMemcpy(d_star, b->star_ps, sizeof(double) * 6 * n, cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_iso, b->iso_ps, sizeof(double) * 6 * n, cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_params, b->params, sizeof(double) * NG_NPARAMS, cudaMemcpyHostToDevice));
    const int threads = 128;
    const int blocks = (n + threads - 1) / threads;
    int slice = NG_FIRST_SLICE < steps_per_launch ? NG_FIRST_SLICE : steps_per_launch;
    int launches = 0;
    for (int s = 0; s < b->n_steps; ) {
        int k = b->n_steps - s < slice ? b->n_steps - s : slice;
        double t0 = ng_now();
        ng_steps<<<blocks, threads>>>(d_star, d_iso, d_state, d_best, d_params, b->dt, s, k, n);
        CK(cudaGetLastError());
        CK(cudaDeviceSynchronize());
        s += k;
        slice = ng_next_slice(k, ng_now() - t0, steps_per_launch);
        launches++;
        if (progress) progress(s, b->n_steps, user);
    }
    fprintf(stderr, "nostos_cuda: %d pairs, %d launches, last slice %d steps\n", n, launches, slice);
    double *best = new double[2 * (size_t)n];
    cudaError_t e = cudaMemcpy(best, d_best, sizeof(double) * 2 * n, cudaMemcpyDeviceToHost);
    if (e == cudaSuccess) {
        for (int i = 0; i < n; i++) { out_min_d[i] = best[2 * i]; out_t_min[i] = best[2 * i + 1]; }
    }
    delete[] best;
    cudaFree(d_star); cudaFree(d_iso); cudaFree(d_state); cudaFree(d_best); cudaFree(d_params);
    if (e != cudaSuccess) { fprintf(stderr, "nostos_cuda: result copy failed: %s\n", cudaGetErrorString(e)); return 1; }
    return 0;
}
