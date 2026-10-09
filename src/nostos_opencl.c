// Nostos OpenCL backend (see nostos_gpu.h). Same device code as the CUDA backend
// (nostos_kernel_body.h, embedded by scripts/embed_kernel.py), double precision required.
// The OpenCL library is loaded at run time (cl_dynload.c, shared with bitboinc's apps); per-vendor
// binaries are built with -DEXPECTED_GPU_VENDOR="amd"/"nvidia" so a task never runs on the wrong GPU.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cl_dynload.h"
#include "nostos_gpu.h"
#include "nostos_kernel_cl.inc"

#ifndef EXPECTED_GPU_VENDOR
#define EXPECTED_GPU_VENDOR ""
#endif

static cl_context g_ctx;
static cl_command_queue g_q;
static cl_program g_prog;
static cl_kernel g_k;
static cl_device_id g_dev;
static size_t g_local = 64;

int ng_cl_init(int platform_id, int device_id) {
    if (!cl_dynload_init()) { fprintf(stderr, "nostos_opencl: could not load an OpenCL library\n"); return 1; }
    cl_uint np = 0;
    if (cl_api.GetPlatformIDs(0, NULL, &np) != CL_SUCCESS || np == 0) { fprintf(stderr, "nostos_opencl: no OpenCL platforms\n"); return 1; }
    cl_platform_id *ps = (cl_platform_id *)malloc(np * sizeof *ps);
    cl_api.GetPlatformIDs(np, ps, NULL);
    // BOINC's --device is an index within the GPU's platform, which need not be platform 0: take the
    // given platform if it matches the vendor and has devices, otherwise the first one that does.
    cl_platform_id plat = 0;
    cl_uint nd = 0;
    if (platform_id >= 0 && (cl_uint)platform_id < np && cl_platform_matches_vendor(ps[platform_id], EXPECTED_GPU_VENDOR)) {
        cl_api.GetDeviceIDs(ps[platform_id], CL_DEVICE_TYPE_GPU, 0, NULL, &nd);
        if (nd > 0) plat = ps[platform_id];
    }
    for (cl_uint i = 0; !plat && i < np; i++) {
        if (!cl_platform_matches_vendor(ps[i], EXPECTED_GPU_VENDOR)) continue;
        cl_uint n = 0;
        if (cl_api.GetDeviceIDs(ps[i], CL_DEVICE_TYPE_GPU, 0, NULL, &n) == CL_SUCCESS && n > 0) { plat = ps[i]; nd = n; }
    }
    free(ps);
    if (!plat) { fprintf(stderr, "nostos_opencl: no \"%s\" OpenCL GPU found\n", EXPECTED_GPU_VENDOR); return 1; }
    cl_device_id *ds = (cl_device_id *)malloc(nd * sizeof *ds);
    cl_api.GetDeviceIDs(plat, CL_DEVICE_TYPE_GPU, nd, ds, NULL);
    if (device_id < 0 || (cl_uint)device_id >= nd) device_id = 0;
    device_id = cl_pick_vendor_device(ds, nd, device_id, EXPECTED_GPU_VENDOR);
    g_dev = ds[device_id];
    free(ds);
    char name[256] = {0};
    cl_api.GetDeviceInfo(g_dev, CL_DEVICE_NAME, sizeof name, name, NULL);
    cl_device_fp_config fp64 = 0;
    cl_api.GetDeviceInfo(g_dev, CL_DEVICE_DOUBLE_FP_CONFIG, sizeof fp64, &fp64, NULL);
    fprintf(stderr, "nostos_opencl: using \"%s\"\n", name);
    if (!fp64) { fprintf(stderr, "nostos_opencl: this GPU has no double precision (fp64), which Nostos needs\n"); return 1; }

    cl_int err;
    g_ctx = cl_api.CreateContext(NULL, 1, &g_dev, NULL, NULL, &err);
    if (!g_ctx) { fprintf(stderr, "nostos_opencl: clCreateContext failed (%d)\n", err); return 1; }
    g_q = cl_api.CreateCommandQueue(g_ctx, g_dev, 0, &err);
    if (!g_q) { fprintf(stderr, "nostos_opencl: clCreateCommandQueue failed (%d)\n", err); return 1; }
    const char *src = kNostosKernelCL;
    size_t len = strlen(src);
    g_prog = cl_api.CreateProgramWithSource(g_ctx, 1, &src, &len, &err);
    err = g_prog ? cl_api.BuildProgram(g_prog, 1, &g_dev, "-cl-std=CL1.2", NULL, NULL) : err;
    if (err != CL_SUCCESS) {
        size_t n = 0;
        cl_api.GetProgramBuildInfo(g_prog, g_dev, CL_PROGRAM_BUILD_LOG, 0, NULL, &n);
        char *log = (char *)malloc(n + 1);
        if (log) {
            cl_api.GetProgramBuildInfo(g_prog, g_dev, CL_PROGRAM_BUILD_LOG, n, log, NULL);
            log[n] = 0;
            fprintf(stderr, "nostos_opencl: kernel build failed (%d):\n%s\n", err, log);
            free(log);
        }
        return 1;
    }
    g_k = cl_api.CreateKernel(g_prog, "ng_steps", &err);
    if (!g_k) { fprintf(stderr, "nostos_opencl: clCreateKernel failed (%d)\n", err); return 1; }
    size_t wg = 0;
    cl_api.GetKernelWorkGroupInfo(g_k, g_dev, CL_KERNEL_WORK_GROUP_SIZE, sizeof wg, &wg, NULL);
    if (wg && wg < g_local) g_local = wg;
    return 0;
}

int ng_cl_run(const NgBatch *b, int steps_per_launch, double *out_min_d, double *out_t_min,
              NgProgress progress, void *user) {
    const int n = b->n_pairs;
    cl_int err;
    int rc = 1;
    cl_mem d_star = cl_api.CreateBuffer(g_ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(double) * 6 * n, (void *)b->star_ps, &err);
    cl_mem d_iso = cl_api.CreateBuffer(g_ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(double) * 6 * n, (void *)b->iso_ps, &err);
    cl_mem d_state = cl_api.CreateBuffer(g_ctx, CL_MEM_READ_WRITE, sizeof(double) * 12 * n, NULL, &err);
    cl_mem d_best = cl_api.CreateBuffer(g_ctx, CL_MEM_READ_WRITE, sizeof(double) * 2 * n, NULL, &err);
    cl_mem d_par = cl_api.CreateBuffer(g_ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(double) * NG_NPARAMS, (void *)b->params, &err);
    double *best = NULL;
    if (!d_star || !d_iso || !d_state || !d_best || !d_par) { fprintf(stderr, "nostos_opencl: buffer allocation failed\n"); goto out; }
    size_t global = ((size_t)n + g_local - 1) / g_local * g_local;
    double dt = b->dt;
    for (int s = 0; s < b->n_steps; s += steps_per_launch) {
        int k = b->n_steps - s < steps_per_launch ? b->n_steps - s : steps_per_launch;
        cl_api.SetKernelArg(g_k, 0, sizeof(cl_mem), &d_star);
        cl_api.SetKernelArg(g_k, 1, sizeof(cl_mem), &d_iso);
        cl_api.SetKernelArg(g_k, 2, sizeof(cl_mem), &d_state);
        cl_api.SetKernelArg(g_k, 3, sizeof(cl_mem), &d_best);
        cl_api.SetKernelArg(g_k, 4, sizeof(cl_mem), &d_par);
        cl_api.SetKernelArg(g_k, 5, sizeof(double), &dt);
        cl_api.SetKernelArg(g_k, 6, sizeof(int), &s);
        cl_api.SetKernelArg(g_k, 7, sizeof(int), &k);
        cl_api.SetKernelArg(g_k, 8, sizeof(int), &n);
        err = cl_api.EnqueueNDRangeKernel(g_q, g_k, 1, NULL, &global, &g_local, 0, NULL, NULL);
        if (err == CL_SUCCESS) err = cl_api.Finish(g_q);
        if (err != CL_SUCCESS) { fprintf(stderr, "nostos_opencl: kernel failed (%d)\n", err); goto out; }
        if (progress) progress(s + k, b->n_steps, user);
    }
    best = (double *)malloc(sizeof(double) * 2 * (size_t)n);
    if (!best || cl_api.EnqueueReadBuffer(g_q, d_best, CL_TRUE, 0, sizeof(double) * 2 * n, best, 0, NULL, NULL) != CL_SUCCESS) {
        fprintf(stderr, "nostos_opencl: result read failed\n"); goto out;
    }
    for (int i = 0; i < n; i++) { out_min_d[i] = best[2 * i]; out_t_min[i] = best[2 * i + 1]; }
    rc = 0;
out:
    free(best);
    if (d_star) cl_api.ReleaseMemObject(d_star);
    if (d_iso) cl_api.ReleaseMemObject(d_iso);
    if (d_state) cl_api.ReleaseMemObject(d_state);
    if (d_best) cl_api.ReleaseMemObject(d_best);
    if (d_par) cl_api.ReleaseMemObject(d_par);
    return rc;
}
