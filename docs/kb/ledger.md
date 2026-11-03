# Pattern ledger

One row per cross-repository pattern ([shared/patterns](shared/patterns)): has this repository been checked for it,
at which commit, and what was found. `kb.py sync` adds an `unchecked` row whenever a new pattern arrives; checking it
means following the pattern's "How to check a repository" section and filling in the row. Status values:
[shared/vocabulary.md](shared/vocabulary.md#ledger). A row that is `affected` or `fixed` names the local issue, and that
issue names the pattern in its `pattern:` field; its Evidence says where, in one clause, and the issue has the rest.
Write rows with `kb.py ledger` (date and commit filled in) rather than by hand.

When code that a row describes changes a lot, check the pattern again and update Checked, Commit and Evidence.

| Pattern | Status | Checked | Commit | Issue | Evidence |
|---|---|---|---|---|---|
| PAT-001 | not-affected | 2026-10-27 | eb7029e |  | CPU: one boinc_time_to_checkpoint() per star, result used once, always paired with boinc_checkpoint_completed() (src/main_nostos_boinc.cpp:143-149). GPU never calls it; its completed() without a begin is clamped at 0 (src/main_nostos_gpu_boinc.cpp:176, see PAT-004) |
| PAT-002 | affected | 2026-10-27 | eb7029e | NOS-001 | start from nostos_checkpoint only, out reopened with "a" and not trimmed (src/main_nostos_boinc.cpp:108-115; GPU src/main_nostos_gpu_boinc.cpp:125-129) |
| PAT-003 | affected | 2026-10-27 | eb7029e | NOS-002 | (a) in place, (b) unchecked, (d) completed() after a failed open: src/main_nostos_boinc.cpp:144-149, src/main_nostos_gpu_boinc.cpp:174-176 |
| PAT-004 | affected | 2026-10-27 | eb7029e | NOS-003 | GPU checkpoint after every batch, without boinc_time_to_checkpoint(), outside a critical section (src/main_nostos_gpu_boinc.cpp:173-176) |
| PAT-005 | affected | 2026-10-27 | eb7029e | NOS-010 | any NULL from boinc_fopen(nostos_checkpoint, "r") means start 0 and out opened with "w": src/main_nostos_boinc.cpp:108-115, src/main_nostos_gpu_boinc.cpp:125-129 |
| PAT-006 | affected | 2026-10-27 | eb7029e | NOS-011 | no fprintf/fflush/ferror/fclose result checked before boinc_finish(0): src/main_nostos_boinc.cpp:136-162, src/main_nostos_gpu_boinc.cpp:169-182 |
| PAT-007 | not-affected | 2026-10-27 | eb7029e |  | CPU has no report at resume, but reports after every star (src/main_nostos_boinc.cpp:152-153) and a star takes about 2.5 s (NOTES 2026-10-08), so the gap is seconds; GPU reports (start + part of batch)/total after the first kernel slice (src/main_nostos_gpu_boinc.cpp:49-53, 153) |
| PAT-008 | n/a | 2026-10-27 | eb7029e |  | no quit_request/abort_request polling, boinc_get_status or direct_process_action in the repository (git grep: no hits); default options (src/main_nostos_boinc.cpp:49-51, src/main_nostos_gpu_boinc.cpp:61-64) |
| PAT-009 | n/a | 2026-10-27 | eb7029e |  | no trickle messages (git grep boinc_send_trickle_up: no hits) |
| PAT-010 | n/a | 2026-10-27 | eb7029e |  | the BOINC worker only computes and reads regular files; no sleep, poll, select or device read in src/ (git grep: no hits); GPU waits are clFinish/cudaDeviceSynchronize |
| PAT-011 | not-affected | 2026-10-27 | eb7029e |  | CUDA: gpu_device_num wins over --device (src/main_nostos_gpu_boinc.cpp:71-72); OpenCL: gpu_opencl_dev_index wins (:75-77); device name logged; this code is under NG_CUDA/NG_OPENCL, defined by every GPU build (build_gpu_server.sh:35, 48, 70; build_gpu_windows.bat:16). An out-of-range index silently becomes 0 (src/nostos_cuda.cu:22, src/nostos_opencl.c:46) |
| PAT-012 | not-affected | 2026-10-27 | eb7029e |  | per-vendor binaries (amd, nvidia; build_gpu_server.sh:49-51, 71-73) take the first platform whose GPU vendor id or vendor name matches the build's vendor, CL_DEVICE_TYPE_GPU for count and list, no CPU fallback: src/nostos_opencl.c:31-48, src/cl_dynload.c:113-165; compiled into every OpenCL build. The build's vendor stands in for gpu_type |
| PAT-013 | affected | 2026-10-27 | eb7029e | NOS-004 | device init and kernel/runtime errors end with boinc_finish(1): src/main_nostos_gpu_boinc.cpp:73, 78, 159 |
| PAT-014 | not-affected | 2026-10-27 | eb7029e |  | CL_DEVICE_DOUBLE_FP_CONFIG checked, device refused with a clear message (src/nostos_opencl.c:52-55); plan classes opencl_ati_fp64/opencl_nvidia_fp64 (dea832d commit message) |
| PAT-015 | not-affected | 2026-10-27 | eb7029e |  | local size 64 clamped to CL_KERNEL_WORK_GROUP_SIZE, which the spec bounds by the device maximum (src/nostos_opencl.c:21, 80-82); global rounded up to a multiple (:98); no __local memory in the kernel |
| PAT-016 | fixed | 2026-10-27 | eb7029e | NOS-012 | fixed 500-step launches hit the Windows TDR (dea832d); since 291f413 slices start at 20 steps and are sized to about 0.25 s: src/nostos_gpu.h:33-67, src/nostos_opencl.c:100-121, src/nostos_cuda.cu:46-58 |
| PAT-017 | n/a | 2026-10-27 | eb7029e | NOS-016 | no Android OpenCL build (scripts/build_boinc_android.sh builds the CPU app only; no Adreno build, dea832d). src/cl_dynload.c lacks planodia 3f6c1b7 (drift recorded in NOS-016) |
| PAT-018 | affected | 2026-10-27 | eb7029e | NOS-013 | 1ecc658 measured different output on the Windows ARM64 and Android builds (upstream clang) |
| PAT-019 | affected | 2026-10-27 | eb7029e | NOS-014 | out opened in text mode in the Windows builds: src/main_nostos_boinc.cpp:115, src/main_nostos_gpu_boinc.cpp:129 |
| PAT-020 | not-affected | 2026-10-27 | eb7029e |  | MinGW builds exist (build_boinc_all.sh:41-43, build_gpu_server.sh:64-76) but src/ creates no threads and has no thread_local (git grep: no hits); only the BOINC API timer thread, which lives until exit |
| PAT-021 | affected | 2026-10-27 | eb7029e | NOS-009, NOS-015 | D:/nostos_data (scripts/download_filter_6d.py:21) and fixed developer paths in scripts/build_boinc_android.sh, build_boinc_linux.sh |
| PAT-022 | affected | 2026-10-27 | eb7029e | NOS-005 | Linux GPU builds in ubuntu:22.04 and nvidia/cuda:12.8.1-devel-ubuntu22.04 (glibc 2.35), glibc linked dynamically: scripts/build_gpu_server.sh:36-38, 52-61 |
| PAT-023 | affected | 2026-10-27 | eb7029e | NOS-006, NOS-008 | scripts/download_filter_6d.py: a file that failed 5 times is marked done (:92-96, :132) |
| PAT-024 | affected | 2026-10-27 | eb7029e | NOS-021 | Linux CUDA 12.8.1: sm_50 to sm_120 plus compute_120 PTX (scripts/build_gpu_server.sh:31) |
| PAT-025 | n/a | 2026-10-27 | eb7029e |  | the BOINC apps start no threads (git grep std::thread, pthread_create, omp, CreateThread, boinc_init_parallel: no hits); only the BOINC API timer thread and GPU driver threads |
| PAT-026 | not-affected | 2026-10-27 | eb7029e |  | the BOINC apps do not sort; qsort only on plain doubles for percentiles in the server tool src/nostos_detail.c:138-140 and in tools/ (ties are equal values, so their order cannot change the output); Python sorted() is stable |
| PAT-027 | affected | 2026-10-27 | eb7029e | NOS-022 | OpenCL: the nine clSetKernelArg results in the launch loop are ignored (src/nostos_opencl.c:105-113) |
| PAT-028 | not-affected | 2026-10-27 | eb7029e |  | guards in shipped code: NG_CUDA/NG_OPENCL with #error when neither is set (src/main_nostos_gpu_boinc.cpp:34-40); EXPECTED_GPU_VENDOR (default empty, which accepts any vendor: src/nostos_opencl.c:12-14) is defined by every OpenCL build (build_gpu_server.sh:51, 73; build_gpu_opencl_wintest.bat:11); the rest are platform guards (_WIN32, __APPLE__) |
| PAT-029 | affected | 2026-10-27 | eb7029e | NOS-023 | GPU builds in nvidia/cuda:12.8.1-devel-ubuntu22.04 and ubuntu:22.04 by tag, no digest (scripts/build_gpu_server.sh:59, 61) |
| PAT-030 | n/a | 2026-10-27 | eb7029e |  | no boinc_temporary_exit() and no deliberate exit without boinc_finish() in src/ (git grep) |
