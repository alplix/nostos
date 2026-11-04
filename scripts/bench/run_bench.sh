#!/bin/bash
# Interleaved speed benchmark of Nostos apps built against the fake BOINC runtime (build_variant.sh).
# usage: run_bench.sh <profile> <rounds> <label>=<binary> [<label>=<binary> ...]
#   profile: a file in scripts/bench/profiles/ without .wu, or a path to a .wu file (its keys: make_bench_wu.py)
#   Runs A B C ... A B C ... (rounds) in fresh run directories, so that load drift hits every binary alike. The input
#   file `in` is made once by make_bench_wu.py and copied into every run directory.
#   env: BENCH_OUT (default $NOSTOS_WORK/bench-runs, NOSTOS_WORK default ~/nostos-work), RUN_PREFIX (e.g.
#        "taskpolicy -b" or "nice -n 19"), NOSTOS_STUB_CHECKPOINT_SECONDS (the stub's checkpoint interval, default 60),
#        GPU_SAMPLE=1 (mean GPU utilisation, sampled every 0.5 s with nvidia-smi, GPU_INDEX picks the card, or ioreg on a
#        Mac; roughly, for the whole GPU). A label may name a wrapper script that sets environment variables (e.g.
#        NOSTOS_BATCH_STARS) and execs the binary.
#   macOS and Linux (GNU time)
# What is compared: the whole `out` file. It holds only the rows `source_id min_pc mean_pc t_myr` and the trailer
# `# done <n_stars>`, no timestamp, timing or host name, so it is the same bytes on every run of one binary. A run whose
# last line is not `# done <n_stars of in>` is marked incomplete.
# Output: one line per run (label, round, wall s, user s, sys s, max RSS MB, exit code, hash of out, mean GPU busy %,
# smallest min_pc in out), then the summary of summarize_bench.py.
set -u
here="$(cd "$(dirname "$0")" && pwd)"
START=$PWD
usage() { awk 'NR>1 && /^#/ {print; next} NR>1 {exit}' "$0"; exit 2; }
[ $# -ge 3 ] || usage
PROFILE=$1; ROUNDS=$2; shift 2
case $ROUNDS in ''|*[!0-9]*) usage ;; esac
WU=$PROFILE; [ -f "$WU" ] || WU=$here/profiles/$PROFILE.wu
[ -f "$WU" ] || { echo "no profile $PROFILE"; exit 2; }
for lb in "$@"; do
  case $lb in *=*) ;; *) echo "'$lb' is not label=binary"; exit 2 ;; esac
  [ -x "${lb#*=}" ] || { echo "${lb#*=} is not an executable file"; exit 2; }
done
STAMP=$(date +%Y%m%d-%H%M%S)
OUT=${BENCH_OUT:-${NOSTOS_WORK:-$HOME/nostos-work}/bench-runs}/$STAMP-$(basename "$WU" .wu)
mkdir -p "$OUT"
RESULTS=$OUT/results.tsv
INPUT=$OUT/in.master
python3 -B "$here/make_bench_wu.py" "$WU" "$INPUT" || { echo "could not make the input of $WU"; exit 2; }
NSTARS=$(awk 'NR>4 && NF==7 {n++} END {print n+0}' "$INPUT")
[ "$NSTARS" -gt 0 ] || { echo "no star lines in $INPUT"; exit 2; }
if [ "$(uname)" = Darwin ]; then
  TIMECMD=(/usr/bin/time -l); RSSDIV=1048576                                   # bytes
else
  TIMECMD=(/usr/bin/time -f "%e real %U user %S sys
%M maximum resident set size"); RSSDIV=1024                                    # GNU time, KiB
fi
# shasum comes with perl, which minimal Linux images lack; sha1sum gives the same hash
HASH=shasum; command -v shasum > /dev/null || HASH=sha1sum
md5f() { md5 -q "$1" 2> /dev/null || md5sum "$1" | cut -d' ' -f1; }
load() { uptime 2> /dev/null | sed 's/.*load/load/' || true; [ -r /proc/loadavg ] && echo "loadavg $(cut -d' ' -f1-3 /proc/loadavg)"; }
gpu_sample() { # one reading in percent
  if command -v nvidia-smi > /dev/null 2>&1; then
    nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader,nounits -i "${GPU_INDEX:-0}" 2> /dev/null | head -1
  elif command -v ioreg > /dev/null 2>&1; then
    ioreg -r -d 1 -c IOAccelerator | grep -o '"Device Utilization %"=[0-9]*' | head -1 | cut -d= -f2
  fi
}
{
  echo "# $(date '+%Y-%m-%d %H:%M:%S') profile $WU rounds $ROUNDS stars $NSTARS host $(sysctl -n machdep.cpu.brand_string 2>/dev/null || uname -m)"
  echo "# in md5 $(md5f "$INPUT"); seed_base $(sed -n 3p "$INPUT"); window, steps and threshold: $(sed -n 4p "$INPUT")"
  echo "# load: $(load | head -1)"
  for lb in "$@"; do echo "# ${lb%%=*} = ${lb#*=} md5 $(md5f "${lb#*=}")"; done
} | tee "$RESULTS"

run_one() { # label binary round
  local label=$1 bin=$2 round=$3 dir=$OUT/$1-$3
  case $bin in /*) ;; *) bin=$START/$bin ;; esac   # the run starts in its own directory
  mkdir -p "$dir"
  cp "$INPUT" "$dir/in"
  local sampler=""
  if [ "${GPU_SAMPLE:-0}" = 1 ]; then
    ( while :; do gpu_sample; sleep 0.5; done ) > "$dir/gpu_util.txt" 2> /dev/null &
    sampler=$!
  fi
  ( cd "$dir" && "${TIMECMD[@]}" ${RUN_PREFIX:-} sh -c 'exec "$@" > stdout.txt 2> stderr.txt' sh "$bin" ) 2> "$dir/time.txt"
  echo $? > "$dir/exit_code"
  [ -n "$sampler" ] && { kill "$sampler"; wait "$sampler"; } 2> /dev/null
  local real user sys rss gpu=- h=none best=-
  real=$(awk '/ real /{print $1}' "$dir/time.txt"); user=$(awk '/ real /{print $3}' "$dir/time.txt"); sys=$(awk '/ real /{print $5}' "$dir/time.txt")
  rss=$(awk -v d=$RSSDIV '/maximum resident set size/{printf "%.0f", $1/d}' "$dir/time.txt")
  [ -s "$dir/gpu_util.txt" ] && gpu=$(awk '$1 ~ /^[0-9.]+$/ {s+=$1; n++} END {if (n) printf "%.0f", s/n}' "$dir/gpu_util.txt")
  if [ -f "$dir/out" ]; then
    # the trailer shows that every star of `in` was read and searched; without it the rows are not a result
    if [ "$(tail -n 1 "$dir/out")" = "# done $NSTARS" ]; then h=$($HASH < "$dir/out" | cut -c1-12); else h=incomplete; fi
    best=$(awk '$1 !~ /^#/ {v = $2 + 0; if (n++ == 0 || v < m) m = v} END {if (n) printf "%.6f", m; else print "-"}' "$dir/out")
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$label" "$round" "$real" "$user" "$sys" "$rss" "$(cat "$dir/exit_code")" "$h" "$gpu" "$best" | tee -a "$RESULTS"
  rm -f "$dir/in" "$dir/nostos_checkpoint"
}
for r in $(seq 1 "$ROUNDS"); do
  for lb in "$@"; do run_one "${lb%%=*}" "${lb#*=}" "$r"; done
done
python3 -B "$here/summarize_bench.py" "$RESULTS" | tee "$OUT/summary.txt"
