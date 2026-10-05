// Minimal 3-vector math for orbit integration. Header-only, inline, no dependencies.
#ifndef NOSTOS_VEC3_H
#define NOSTOS_VEC3_H

#include <math.h>

typedef struct { double x, y, z; } Vec3;

static inline Vec3 v3_add(Vec3 a, Vec3 b) { return (Vec3){a.x+b.x, a.y+b.y, a.z+b.z}; }
static inline Vec3 v3_sub(Vec3 a, Vec3 b) { return (Vec3){a.x-b.x, a.y-b.y, a.z-b.z}; }
static inline Vec3 v3_scale(Vec3 a, double s) { return (Vec3){a.x*s, a.y*s, a.z*s}; }
static inline double v3_dot(Vec3 a, Vec3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static inline double v3_norm(Vec3 a) { return sqrt(v3_dot(a, a)); }

static inline Vec3 v3_cross(Vec3 a, Vec3 b) {
    return (Vec3){
        a.y*b.z - a.z*b.y,
        a.z*b.x - a.x*b.z,
        a.x*b.y - a.y*b.x
    };
}

#endif
