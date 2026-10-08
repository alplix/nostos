#!/bin/bash
# Static Linux x86_64 build of the Nostos CPU BOINC app (src/main_nostos_boinc.cpp), run inside WSL
# (Ubuntu) on the dev PC. Links against the same locally built BOINC 8.2.8 libraries the sibling
# bitboinc-app project uses for its local Linux builds (see bitboinc-app/build_potamos_gpu_linux.sh):
# sources /root/boincsrc/boinc, out-of-tree build /root/boinc-build-x86 (config.h + .libs/*.a).
# The very first build (2026-10-06) was linked on bitboinc-app's server against
# /opt/bitboinc-src/boinc; that tree is not used here (no SSH from this machine).
#
# Usage (from Windows Git Bash):
#   MSYS_NO_PATHCONV=1 wsl -d Ubuntu -u root -- bash /mnt/c/Users/Alp/Documents/nostos/scripts/build_boinc_linux.sh
# Output: build/linux/nostos_cpu_linux_x86_64 (git-ignored). Runs at nice 19.
set -euo pipefail
REPO=/mnt/c/Users/Alp/Documents/nostos
BSRC=${BOINC_SRC:-/root/boincsrc/boinc}
BBLD=${BOINC_BUILD:-/root/boinc-build-x86}
OUT=$REPO/build/linux
mkdir -p "$OUT"
cd "$REPO"

CSRC="src/coords.c src/orbital_elements.c src/potential.c src/encounter.c src/montecarlo.c src/mc_search.c src/nostos_io.c"
OBJ=""
TMP=$(mktemp -d)
for c in $CSRC; do
    o="$TMP/$(basename "${c%.c}").o"
    nice -n 19 gcc -O2 -std=c99 -Wall -Wextra -c "$c" -o "$o"
    OBJ="$OBJ $o"
done
nice -n 19 g++ -O2 -Wall -Wextra -static \
    -I src -I "$BSRC/api" -I "$BSRC/lib" -I "$BBLD" \
    src/main_nostos_boinc.cpp $OBJ \
    "$BBLD/api/.libs/libboinc_api.a" "$BBLD/lib/.libs/libboinc.a" \
    -lpthread -lm -o "$OUT/nostos_cpu_linux_x86_64"
rm -rf "$TMP"
file "$OUT/nostos_cpu_linux_x86_64"
echo "glibc-versioned dynamic symbols: $(objdump -T "$OUT/nostos_cpu_linux_x86_64" 2>/dev/null | grep -c GLIBC_ || true)"
ls -la "$OUT/nostos_cpu_linux_x86_64"
