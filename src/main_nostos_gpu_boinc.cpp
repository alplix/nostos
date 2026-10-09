// Nostos GPU BOINC app (CUDA or OpenCL, chosen at compile time: -DNG_CUDA or -DNG_OPENCL).
//
// Same work-unit input and output format as the CPU app (src/main_nostos_boinc.cpp): an ISO header
// and one star per line in, "source_id min_pc mean_pc t_myr" for every star under the threshold out,
// then "# done <n_stars>". The physics runs on the GPU (nostos_gpu.h); everything that decides the
// numbers on the host side -- ISO start state, star phase space, clone draws, per-star reduction --
// is the CPU app's own code in the CPU app's order, so GPU and CPU results of one work unit agree.
//
// Stars are processed in batches (NG_BATCH_STARS, default 2000; env NOSTOS_BATCH_STARS); the
// checkpoint is the index of the next star, written after each finished batch (a batch is a minute
// or two on a mid-range GPU). Each kernel launch runs a slice of the 100,000 steps
// sized to take ~0.25 s (nostos_gpu.h; NOSTOS_STEPS_PER_LAUNCH caps it, default 5000).
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <vector>

#include "boinc_api.h"
#include "filesys.h"

extern "C" {
#include "coords.h"
#include "orbital_elements.h"
#include "potential.h"
#include "montecarlo.h"
#include "mc_search.h"
#include "units.h"
#include "nostos_io.h"
}
#include "nostos_gpu.h"

#if defined(NG_CUDA)
#define NG_NAME "nostos_cuda"
#elif defined(NG_OPENCL)
#define NG_NAME "nostos_opencl"
#else
#error "build with -DNG_CUDA or -DNG_OPENCL"
#endif

static int env_int(const char *name, int def) {
    const char *v = getenv(name);
    int x = v ? atoi(v) : 0;
    return x > 0 ? x : def;
}

struct Progress { size_t stars_before, n_stars, batch_stars; };
static void on_progress(int steps_done, int n_steps, void *user) {
    Progress *p = (Progress *)user;
    double f = (p->stars_before + (double)p->batch_stars * steps_done / n_steps) / (double)p->n_stars;
    boinc_fraction_done(f);
}

static int arg_int(int argc, char **argv, const char *name, int def) {
    for (int i = 1; i + 1 < argc; i++) if (!strcmp(argv[i], name)) return atoi(argv[i + 1]);
    return def;
}

int main(int argc, char **argv) {
    BOINC_OPTIONS options;
    boinc_options_defaults(options);
    options.normal_thread_priority = true;  // keep the GPU fed
    int retval = boinc_init_options(&options);
    if (retval) { fprintf(stderr, "boinc_init_options failed: %d\n", retval); exit(retval); }

    // device: app_init_data first (BOINC's own assignment), --device as a fallback
    APP_INIT_DATA aid;
    boinc_get_init_data(aid);
#if defined(NG_CUDA)
    int device = arg_int(argc, argv, "--device", 0);
    if (aid.gpu_device_num >= 0) device = aid.gpu_device_num;
    if (ng_cuda_init(device)) { boinc_finish(1); return 1; }
#else
    int device = arg_int(argc, argv, "--device", 0);
    int platform = arg_int(argc, argv, "--cl-platform", 0);
    if (aid.gpu_opencl_dev_index >= 0) device = aid.gpu_opencl_dev_index;
    if (ng_cl_init(platform, device)) { boinc_finish(1); return 1; }
#endif

    char input_path[512], output_path[512];
    if (boinc_resolve_filename("in", input_path, sizeof(input_path)) ||
        boinc_resolve_filename("out", output_path, sizeof(output_path))) {
        fprintf(stderr, "can't resolve in/out\n"); boinc_finish(1); return 1;
    }
    OrbitalElements el;
    double vel_cov[3][3];
    uint64_t seed_base;
    double total_time_myr, threshold_pc;
    int n_steps;
    StarRecord *stars = NULL;
    size_t n_stars = 0;
    FILE *inf = boinc_fopen(input_path, "r");
    bool parsed = inf && parse_input_stream(inf, &el, vel_cov, &seed_base, &total_time_myr, &n_steps,
                                             &threshold_pc, &stars, &n_stars);
    if (inf) fclose(inf);
    if (!parsed) { fprintf(stderr, "can't open or parse %s\n", input_path); boinc_finish(1); return 1; }

    // ISO start: exactly as the CPU app (inbound asymptote, t = 0)
    PhaseSpace iso_ps;
    if (!orbital_elements_inbound_asymptote_gm(el, NOSTOS_GM_SSB, 0.0, &solar_params_default, &iso_ps)) {
        fprintf(stderr, "invalid orbital elements\n"); free(stars); boinc_finish(1); return 1;
    }
    // Cholesky of the clone covariance: the same for every star (mc_search_stage1)
    double cov6[MC_MAX_DIM][MC_MAX_DIM] = {{0}};
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) cov6[i][j] = vel_cov[i][j];
    double L[MC_MAX_DIM][MC_MAX_DIM];
    const bool cov_ok = mc_cholesky(cov6, 3, L) == 0;
    if (!cov_ok) fprintf(stderr, "covariance not positive definite: no star can be searched (same as the CPU app)\n");
    const double mean6[MC_MAX_DIM] = {iso_ps.vel.x, iso_ps.vel.y, iso_ps.vel.z, 0, 0, 0};

    NgBatch b;
    const MilkyWayParams &mw = milkyway_potential_default;
    const double prm[NG_NPARAMS] = {mw.bulge_GM, mw.bulge_b, mw.disk_GM, mw.disk_a, mw.disk_b,
                                    mw.halo_GM, mw.halo_a, mw.halo_cutoff, mw.halo_gamma};
    memcpy(b.params, prm, sizeof(prm));
    b.dt = total_time_myr / n_steps;
    b.n_steps = n_steps;
    const double threshold_kpc = threshold_pc / 1000.0;
    const int C = MC_N_ISO_CLONES;
    const size_t batch = (size_t)env_int("NOSTOS_BATCH_STARS", 2000);
    const int steps_per_launch = env_int("NOSTOS_STEPS_PER_LAUNCH", 5000);  // upper bound; actual slices adapt (nostos_gpu.h)

    size_t start = 0;
    FILE *cp = boinc_fopen("nostos_checkpoint", "r");
    if (cp) { if (fscanf(cp, "%zu", &start) != 1) start = 0; fclose(cp); }
    if (start > n_stars) start = 0;
    if (start) fprintf(stderr, "resuming from checkpoint: star %zu of %zu\n", start, n_stars);
    FILE *outf = boinc_fopen(output_path, start > 0 ? "a" : "w");
    if (!outf) { fprintf(stderr, "can't open %s\n", output_path); free(stars); boinc_finish(1); return 1; }

    for (size_t s0 = cov_ok ? start : n_stars; s0 < n_stars; s0 += batch) {
        size_t ns = n_stars - s0 < batch ? n_stars - s0 : batch;
        std::vector<double> star_ps(ns * C * 6), iso_v(ns * C * 6), min_d(ns * C), t_min(ns * C);
        for (size_t k = 0; k < ns; k++) {
            const StarRecord &st = stars[s0 + k];
            PhaseSpace sp;
            astrometry_to_galactocentric(st.astrometry, &solar_params_default, &sp);
            McRng rng;
            mc_rng_seed(&rng, combine_seed(seed_base, st.source_id));
            for (int c = 0; c < C; c++) {
                double drawn[MC_MAX_DIM];
                mc_sample_multivariate_gaussian(&rng, mean6, L, 3, drawn);
                double *o = &star_ps[(k * C + c) * 6];
                o[0] = sp.pos.x; o[1] = sp.pos.y; o[2] = sp.pos.z; o[3] = sp.vel.x; o[4] = sp.vel.y; o[5] = sp.vel.z;
                double *q = &iso_v[(k * C + c) * 6];
                q[0] = iso_ps.pos.x; q[1] = iso_ps.pos.y; q[2] = iso_ps.pos.z; q[3] = drawn[0]; q[4] = drawn[1]; q[5] = drawn[2];
            }
        }
        b.n_pairs = (int)(ns * C);
        b.star_ps = star_ps.data();
        b.iso_ps = iso_v.data();
        Progress pr = {s0, n_stars, ns};
#if defined(NG_CUDA)
        int rc = ng_cuda_run(&b, steps_per_launch, min_d.data(), t_min.data(), on_progress, &pr);
#else
        int rc = ng_cl_run(&b, steps_per_launch, min_d.data(), t_min.data(), on_progress, &pr);
#endif
        if (rc) { fprintf(stderr, NG_NAME ": GPU run failed\n"); fclose(outf); free(stars); boinc_finish(1); return 1; }
        // per-star reduction in clone order, exactly as mc_search_stage1
        for (size_t k = 0; k < ns; k++) {
            double best = INFINITY, sum = 0.0, t_best = 0.0;
            for (int c = 0; c < C; c++) {
                double d = min_d[k * C + c];
                sum += d;
                if (d < best) { best = d; t_best = t_min[k * C + c]; }
            }
            if (best < threshold_kpc) {
                fprintf(outf, "%llu %.6f %.6f %.6f\n", (unsigned long long)stars[s0 + k].source_id,
                        best * 1000.0, (sum / C) * 1000.0, t_best);
            }
        }
        fflush(outf);
        FILE *cpw = boinc_fopen("nostos_checkpoint", "w");
        if (cpw) { fprintf(cpw, "%zu\n", s0 + ns); fclose(cpw); }
        boinc_checkpoint_completed();
        boinc_fraction_done((double)(s0 + ns) / n_stars);
    }
    fprintf(outf, "# done %zu\n", n_stars);
    fclose(outf);
    free(stars);
    boinc_finish(0);
    return 0;
}
