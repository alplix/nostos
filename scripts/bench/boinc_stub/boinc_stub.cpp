// Fake BOINC runtime for running the Nostos wrappers as plain processes in scripts/bench/ (see boinc_api.h).
// Files are plain files in the working directory, so boinc_resolve_filename("in") is "./in".
//
// What differs from the real library, and so from a volunteer's run:
//  * boinc_finish() exits at once; the real one sleeps 2 s first (docs/kb/shared/boinc-api.md). A run here is 2 s
//    shorter than a task, which only matters for the short profiles.
//  * No timer thread, no heartbeat, no suspend, quit or abort, no priority change (the client lowers the priority of an
//    app; scripts/bench/run_bench.sh has RUN_PREFIX for that).
//  * boinc_fraction_done() writes nothing (the real one updates shared memory, about once a second at most).
//  * boinc_time_to_checkpoint() follows the real contract (1 at most once per period, re-armed by
//    boinc_checkpoint_completed()) but enters no critical section. The period is NOSTOS_STUB_CHECKPOINT_SECONDS
//    (default 60, the client's default; 0 = every call). The CPU wrapper writes nostos_checkpoint when it gets a 1; the
//    GPU wrapper writes it after every batch whatever this says (NOS-003).
#include "boinc_api.h"
#include "filesys.h"
#include <cstdio>
#include <cstdlib>
#include <time.h>

static double stub_period = 60.0;
static double stub_last = 0.0;

static double stub_now() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + ts.tv_nsec * 1e-9;
}

int boinc_init_options(BOINC_OPTIONS*) {
    const char* p = getenv("NOSTOS_STUB_CHECKPOINT_SECONDS");
    if (p) stub_period = atof(p);
    stub_last = stub_now();
    return 0;
}

int boinc_finish(int status) {
    fprintf(stderr, "boinc_stub: boinc_finish(%d)\n", status);
    fflush(stderr);
    exit(status);
}

int boinc_resolve_filename(const char* virtual_name, char* physical_name, int len) {
    snprintf(physical_name, len, "%s", virtual_name);
    return 0;
}

int boinc_get_init_data(APP_INIT_DATA& d) {
    d = APP_INIT_DATA();
    return 0;
}

int boinc_fraction_done(double) { return 0; }

int boinc_time_to_checkpoint() {
    return stub_now() - stub_last >= stub_period ? 1 : 0;
}

int boinc_checkpoint_completed() {
    stub_last = stub_now();
    return 0;
}

FILE* boinc_fopen(const char* path, const char* mode) { return fopen(path, mode); }
