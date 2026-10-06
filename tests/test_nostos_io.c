// Validates nostos_io.c's work-unit input parsing (main_nostos_boinc.cpp's own I/O logic,
// pulled out specifically so it could be tested without needing BOINC headers/libs locally --
// see nostos_io.h's top comment).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../src/nostos_io.h"

static int g_failures = 0;

static void check(int ok, const char *what) {
    if (!ok) { printf("  FAIL: %s\n", what); g_failures++; }
    else { printf("  ok:   %s\n", what); }
}

static FILE *open_string_as_stream(const char *s) {
    // tmpfile() + write-then-rewind: portable way to get a real FILE* backed by this exact
    // content, so parse_input_stream() is tested through the identical fscanf-based code path
    // main_nostos_boinc.cpp actually uses (not a hand-rolled string-parsing stand-in).
    FILE *f = tmpfile();
    fputs(s, f);
    rewind(f);
    return f;
}

static void test_parses_valid_input(void) {
    printf("test_parses_valid_input:\n");
    const char *input =
        "1.2 0.26 122.7 24.6 241.8\n"
        "1e-4 0.0 0.0 1e-4 0.0 1e-4\n"
        "12345\n"
        "-2.0 20000 50.0\n"
        "2525688198820543360 12.056 -5.136 40.2 182.16 -120.16 26.88\n"
        "3107000885484340224 101.095 -1.853 24.12 -5.33 -70.91 33.12\n";
    FILE *f = open_string_as_stream(input);

    OrbitalElements el;
    double vel_cov[3][3];
    uint64_t seed_base;
    double total_time_myr, distance_threshold_pc;
    int n_steps;
    StarRecord *stars = NULL;
    size_t n_stars = 0;
    bool ok = parse_input_stream(f, &el, vel_cov, &seed_base, &total_time_myr, &n_steps,
                                  &distance_threshold_pc, &stars, &n_stars);
    fclose(f);

    check(ok, "parses a well-formed input stream");
    check(fabs(el.e - 1.2) < 1e-12 && fabs(el.q_au - 0.26) < 1e-12, "orbital elements read correctly");
    check(fabs(vel_cov[0][0] - 1e-4) < 1e-12 && fabs(vel_cov[1][1] - 1e-4) < 1e-12 &&
          fabs(vel_cov[0][1] - 0.0) < 1e-12, "velocity covariance read correctly (diagonal and off-diagonal)");
    check(vel_cov[1][0] == vel_cov[0][1] && vel_cov[2][0] == vel_cov[0][2] && vel_cov[2][1] == vel_cov[1][2],
          "velocity covariance is correctly mirrored into a full symmetric matrix from the upper-triangle input");
    check(seed_base == 12345ULL, "seed_base read correctly");
    check(total_time_myr == -2.0 && n_steps == 20000 && distance_threshold_pc == 50.0,
          "search parameters read correctly");
    check(n_stars == 2, "reads both star records");
    if (n_stars == 2) {
        check(stars[0].source_id == 2525688198820543360ULL, "first star's source_id (a large uint64) read correctly");
        check(fabs(stars[0].astrometry.ra_deg - 12.056) < 1e-9, "first star's astrometry read correctly");
        check(stars[1].source_id == 3107000885484340224ULL, "second star's source_id read correctly");
    }
    free(stars);
}

static void test_rejects_malformed_input(void) {
    printf("test_rejects_malformed_input:\n");
    const char *input = "not_a_number 0.26 122.7 24.6 241.8\n";
    FILE *f = open_string_as_stream(input);

    OrbitalElements el;
    double vel_cov[3][3];
    uint64_t seed_base;
    double total_time_myr, distance_threshold_pc;
    int n_steps;
    StarRecord *stars = NULL;
    size_t n_stars = 0;
    bool ok = parse_input_stream(f, &el, vel_cov, &seed_base, &total_time_myr, &n_steps,
                                  &distance_threshold_pc, &stars, &n_stars);
    fclose(f);
    check(!ok, "rejects a file that can't even parse the orbital elements line");
}

static void test_zero_stars_is_valid(void) {
    printf("test_zero_stars_is_valid:\n");
    // A work unit with header lines but no star data at all -- should parse cleanly with
    // n_stars=0, not be treated as an error (an empty chunk is a legitimate, if unusual, input).
    const char *input =
        "1.2 0.26 122.7 24.6 241.8\n"
        "1e-4 0.0 0.0 1e-4 0.0 1e-4\n"
        "12345\n"
        "-2.0 20000 50.0\n";
    FILE *f = open_string_as_stream(input);

    OrbitalElements el;
    double vel_cov[3][3];
    uint64_t seed_base;
    double total_time_myr, distance_threshold_pc;
    int n_steps;
    StarRecord *stars = NULL;
    size_t n_stars = 0;
    bool ok = parse_input_stream(f, &el, vel_cov, &seed_base, &total_time_myr, &n_steps,
                                  &distance_threshold_pc, &stars, &n_stars);
    fclose(f);
    check(ok, "parses successfully even with zero star records");
    check(n_stars == 0, "n_stars is 0, not an error");
    free(stars);
}

static void test_many_stars_grows_correctly(void) {
    printf("test_many_stars_grows_correctly:\n");
    // 200 star records -- well past the initial 64-element capacity in parse_input_stream()'s
    // malloc/realloc growth, exercising the realloc path (and that it doesn't corrupt or drop
    // earlier entries when it grows).
    char buf[1 << 16];
    int pos = snprintf(buf, sizeof(buf), "1.2 0.26 122.7 24.6 241.8\n1e-4 0 0 1e-4 0 1e-4\n1\n-2.0 1000 50.0\n");
    for (int i = 0; i < 200; i++) {
        pos += snprintf(buf + pos, sizeof(buf) - pos, "%d 10.0 -5.0 40.0 100.0 -50.0 20.0\n", 1000 + i);
    }
    FILE *f = open_string_as_stream(buf);

    OrbitalElements el;
    double vel_cov[3][3];
    uint64_t seed_base;
    double total_time_myr, distance_threshold_pc;
    int n_steps;
    StarRecord *stars = NULL;
    size_t n_stars = 0;
    bool ok = parse_input_stream(f, &el, vel_cov, &seed_base, &total_time_myr, &n_steps,
                                  &distance_threshold_pc, &stars, &n_stars);
    fclose(f);
    check(ok, "parses successfully with 200 star records (past the initial realloc capacity)");
    check(n_stars == 200, "all 200 stars were read, none dropped during growth");
    if (n_stars == 200) {
        check(stars[0].source_id == 1000ULL, "first star survived the later reallocations");
        check(stars[199].source_id == 1199ULL, "last star (the one that triggered the final growth) is correct");
    }
    free(stars);
}

static void test_combine_seed_properties(void) {
    printf("test_combine_seed_properties:\n");
    uint64_t s1 = combine_seed(42, 1000);
    uint64_t s2 = combine_seed(42, 1000);
    uint64_t s3 = combine_seed(42, 1001);
    uint64_t s4 = combine_seed(43, 1000);
    check(s1 == s2, "same (seed_base, source_id) always produces the same combined seed (determinism)");
    check(s1 != s3, "a different source_id produces a different combined seed");
    check(s1 != s4, "a different seed_base produces a different combined seed");
}

int main(void) {
    test_parses_valid_input();
    test_rejects_malformed_input();
    test_zero_stars_is_valid();
    test_many_stars_grows_correctly();
    test_combine_seed_properties();

    if (g_failures) {
        printf("\n%d check(s) FAILED\n", g_failures);
        return 1;
    }
    printf("\nAll checks passed.\n");
    return 0;
}
