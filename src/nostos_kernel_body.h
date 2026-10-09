// Shared device code for the Nostos GPU backends: included by nostos_cuda.cu (as CUDA) and embedded
// as source text by nostos_opencl.c (as OpenCL C). Written in the common subset of both languages:
// NG_GLOBAL / NG_DEVICE / NG_KERNEL / NG_TID are defined by the includer.
//
// Literal port of src/potential.c milkyway_accel() and src/integrator.h leapfrog_step(), operation by
// operation (no algebraic rewrites, FP contraction off), so results stay within rounding of the CPU.

NG_DEVICE void ng_accel(const double *p, double x, double y, double z, double *ax, double *ay, double *az) {
    double R2 = x * x + y * y;
    double R = sqrt(R2);
    double r = sqrt(x * x + y * y + z * z);
    double tx = 0.0, ty = 0.0, tz = 0.0;
    {   // Plummer bulge
        double denom = pow(r * r + p[1] * p[1], 1.5);
        double f = -p[0] / (denom > 0 ? denom : 1e-12);
        tx = tx + x * f; ty = ty + y * f; tz = tz + z * f;
    }
    {   // Miyamoto-Nagai disk
        double zeta = sqrt(z * z + p[4] * p[4]);
        double azeta = p[3] + zeta;
        double D2 = R2 + azeta * azeta;
        double D = sqrt(D2);
        double D3 = D2 * D;
        double aR = -p[2] * R / D3;
        double azz = -p[2] * z * azeta / (zeta > 0 ? zeta : 1e-12) / D3;
        if (R > 1e-12) {
            tx += aR * x / R;
            ty += aR * y / R;
        }
        tz += azz;
    }
    {   // truncated Allen & Santillan halo
        double r_safe = (r > 1e-12 ? r : 1e-12);
        double r_eff = (r_safe < p[7]) ? r_safe : p[7];
        double xx = r_eff / p[6];
        double gm = p[5] * pow(xx, p[8]) / (1.0 + pow(xx, p[8] - 1.0));
        double f = -gm / (r_safe * r_safe * r_safe);
        tx = tx + x * f; ty = ty + y * f; tz = tz + z * f;
    }
    *ax = tx; *ay = ty; *az = tz;
}

// One leapfrog step of one body: s = pos.xyz, vel.xyz.
NG_DEVICE void ng_leapfrog(const double *p, double *s, double dt) {
    double ax, ay, az;
    ng_accel(p, s[0], s[1], s[2], &ax, &ay, &az);
    double h = dt * 0.5;
    double vx = s[3] + ax * h, vy = s[4] + ay * h, vz = s[5] + az * h;
    s[0] = s[0] + vx * dt; s[1] = s[1] + vy * dt; s[2] = s[2] + vz * dt;
    ng_accel(p, s[0], s[1], s[2], &ax, &ay, &az);
    s[3] = vx + ax * h; s[4] = vy + ay * h; s[5] = vz + az * h;
}

// state per pair: 12 doubles (star pos/vel, ISO pos/vel); best per pair: min distance, time.
// Runs steps [step0, step0 + nsteps) of every pair; the first call (step0 == 0) also initialises
// state from the initial phase spaces and the minimum from the t=0 separation, exactly like
// find_closest_approach().
NG_KERNEL void ng_steps(NG_GLOBAL const double *star_ps, NG_GLOBAL const double *iso_ps,
                        NG_GLOBAL double *state, NG_GLOBAL double *best,
                        NG_GLOBAL const double *params, double dt, int step0, int nsteps, int n_pairs) {
    int i = NG_TID;
    if (i >= n_pairs) return;
    double p[9];
    for (int k = 0; k < 9; k++) p[k] = params[k];
    double st[6], is[6];
    double min_d, t_min;
    if (step0 == 0) {
        for (int k = 0; k < 6; k++) { st[k] = star_ps[i * 6 + k]; is[k] = iso_ps[i * 6 + k]; }
        double sx = st[0] - is[0], sy = st[1] - is[1], sz = st[2] - is[2];
        min_d = sqrt(sx * sx + sy * sy + sz * sz);
        t_min = 0.0;
    } else {
        for (int k = 0; k < 6; k++) { st[k] = state[i * 12 + k]; is[k] = state[i * 12 + 6 + k]; }
        min_d = best[i * 2];
        t_min = best[i * 2 + 1];
    }
    for (int j = step0; j < step0 + nsteps; j++) {
        ng_leapfrog(p, st, dt);
        ng_leapfrog(p, is, dt);
        double sx = st[0] - is[0], sy = st[1] - is[1], sz = st[2] - is[2];
        double d = sqrt(sx * sx + sy * sy + sz * sz);
        if (d < min_d) {
            min_d = d;
            t_min = (j + 1) * dt;
        }
    }
    for (int k = 0; k < 6; k++) { state[i * 12 + k] = st[k]; state[i * 12 + 6 + k] = is[k]; }
    best[i * 2] = min_d;
    best[i * 2 + 1] = t_min;
}
