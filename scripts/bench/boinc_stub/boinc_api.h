#ifndef NOSTOS_BENCH_BOINC_API_H
#define NOSTOS_BENCH_BOINC_API_H
// Minimal stand-in for BOINC's api/boinc_api.h: enough to build the two Nostos BOINC wrappers
// (src/main_nostos_boinc.cpp, src/main_nostos_gpu_boinc.cpp) without a BOINC checkout and run them as plain processes,
// for scripts/bench/. Only what the wrappers call is declared; boinc_stub.cpp says what each function does here. This is
// a speed-benchmark runtime, not a test harness: it has no kill and resume support (docs/BENCHMARK.md).
//
// Written without C++11 features (no default member initialisers): the CPU wrapper is built with the compiler's default
// language level, as scripts/build_boinc_all.sh builds it, and Apple clang's default is C++98.

struct BOINC_OPTIONS {
    int main_program;
    int check_heartbeat;
    int handle_process_control;
    int send_status_msgs;
    int direct_process_action;
    int multi_thread;
    int multi_process;
    int normal_thread_priority;
    BOINC_OPTIONS()
        : main_program(0), check_heartbeat(0), handle_process_control(0), send_status_msgs(0),
          direct_process_action(0), multi_thread(0), multi_process(0), normal_thread_priority(0) {}
};

// -1 = the client assigned no device, so the GPU wrapper falls back to --device (default 0)
struct APP_INIT_DATA {
    int gpu_device_num;
    int gpu_opencl_dev_index;
    double checkpoint_period;
    APP_INIT_DATA() : gpu_device_num(-1), gpu_opencl_dev_index(-1), checkpoint_period(60) {}
};

inline void boinc_options_defaults(BOINC_OPTIONS& o) { o = BOINC_OPTIONS(); }
int boinc_init_options(BOINC_OPTIONS*);
int boinc_finish(int status);
int boinc_resolve_filename(const char* virtual_name, char* physical_name, int len);
int boinc_get_init_data(APP_INIT_DATA&);
int boinc_fraction_done(double);
int boinc_time_to_checkpoint();
int boinc_checkpoint_completed();

#endif
