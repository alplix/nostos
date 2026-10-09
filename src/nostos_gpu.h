// Nostos GPU backends (CUDA: nostos_cuda.cu, OpenCL: nostos_opencl.c) -- one interface, used by
// main_nostos_gpu_boinc.cpp.
//
// The GPU runs exactly the CPU app's stage-1 search (mc_search_stage1 -> find_closest_approach ->
// leapfrog_step with the Allen & Santillan potential, src/potential.c), in double precision, one
// thread per (star, ISO clone) pair. The clone velocities are drawn on the host with the CPU app's
// own RNG and Cholesky (src/montecarlo.c), in the CPU app's order, and the per-star reduction
// (min over clones, mean = sum/100 in clone order) is done on the host the same way -- so a GPU
// result agrees with a CPU result of the same work unit to rounding level, which quorum-2
// validation (server/nostos_compare.sh) relies on.
#ifndef NOSTOS_GPU_H
#define NOSTOS_GPU_H

#ifdef __cplusplus
extern "C" {
#endif

// Potential parameters, in the order the kernels expect (MilkyWayParams, src/potential.h).
enum { NG_BULGE_GM, NG_BULGE_B, NG_DISK_GM, NG_DISK_A, NG_DISK_B, NG_HALO_GM, NG_HALO_A, NG_HALO_CUTOFF, NG_HALO_GAMMA, NG_NPARAMS };

typedef struct {
    int n_pairs;               // stars x clones
    const double *star_ps;     // per pair: star pos.xyz, vel.xyz (6 doubles)
    const double *iso_ps;      // per pair: ISO clone pos.xyz, vel.xyz (6 doubles)
    double params[NG_NPARAMS];
    double dt;                 // Myr (negative: backwards)
    int n_steps;
} NgBatch;

// progress(steps_done, n_steps, user) is called after every slice of steps.
typedef void (*NgProgress)(int steps_done, int n_steps, void *user);

// CUDA backend. Returns 0 on success.
int ng_cuda_init(int device);
int ng_cuda_run(const NgBatch *b, int steps_per_launch, double *out_min_d, double *out_t_min,
                NgProgress progress, void *user);

// OpenCL backend (needs a device with double precision). Returns 0 on success.
int ng_cl_init(int platform, int device);
int ng_cl_run(const NgBatch *b, int steps_per_launch, double *out_min_d, double *out_t_min,
              NgProgress progress, void *user);

#ifdef __cplusplus
}
#endif
#endif
