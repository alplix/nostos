#!/bin/bash
# One worker of a split full-catalogue scan (src/gaia_scan.c --byte-range). Start N copies (N = the
# thread budget); they claim the NCHUNK byte ranges of the CSV dynamically (atomic mkdir), so an
# uneven candidate density along the file (the CSV is in source_id = sky order) cannot leave one
# worker idle while the other still has half the file. Each chunk writes <out_prefix>_cNN.txt/.log.
# Launch it at Windows Idle priority (child processes inherit the priority class), e.g. from
# PowerShell:  $p = Start-Process bash -ArgumentList "scripts/scan_worker.sh", ... -PassThru;
#              $p.PriorityClass = 'Idle'
# Usage: scan_worker.sh <gaia_scan_exe> <csv> <out_prefix> <nchunk> [extra gaia_scan args...]
set -u
EXE=$1; CSV=$2; OUT=$3; NCHUNK=$4; shift 4
SIZE=$(stat -c %s "$CSV")
CLAIMS="${OUT}_claims"
mkdir -p "$CLAIMS"
for ((i = 0; i < NCHUNK; i++)); do
    mkdir "$CLAIMS/$i" 2>/dev/null || continue   # someone else has it
    START=$(( SIZE * i / NCHUNK ))
    END=$(( SIZE * (i + 1) / NCHUNK ))
    n=$(printf "%02d" "$i")
    "$EXE" "$CSV" "$@" --byte-range=$START:$END > "${OUT}_c$n.txt" 2> "${OUT}_c$n.log"
    echo "chunk $i exit $? $(date +%H:%M:%S)" >> "${OUT}_workers.log"
done
