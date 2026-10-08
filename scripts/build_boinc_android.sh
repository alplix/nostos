#!/bin/bash
# Nostos CPU BOINC app for Android aarch64 (WSL "Ubuntu", as root):
#   MSYS_NO_PATHCONV=1 wsl -d Ubuntu -u root -- bash /mnt/c/Users/Alp/Documents/nostos/scripts/build_boinc_android.sh
# Same NDK r27c / API 21 BOINC libs as bitboinc's Android builds (/root/boinc-build-ndk). libc++ is
# linked statically: Android doesn't ship libc++_shared.so to apps (bitboinc forum thread 74).
set -e
V=${NOSTOS_VERSION:-1}
SRC=/mnt/c/Users/Alp/Documents/nostos
OUT=$SRC/build/android
NDK=/root/android/android-ndk-r27c/toolchains/llvm/prebuilt/linux-x86_64/bin
CC=$NDK/aarch64-linux-android21-clang
CXX=$NDK/aarch64-linux-android21-clang++
BOINC=/root/boinc-build-ndk
BOINCSRC=/root/wb063/boinc-auto
W=/tmp/nostos_android; rm -rf $W; mkdir -p $W $OUT
for c in coords orbital_elements potential encounter montecarlo mc_search nostos_io; do
    $CC -O2 -std=c99 -I $SRC/src -c $SRC/src/$c.c -o $W/$c.o
done
$CXX -O2 -I $SRC/src -I $BOINCSRC/api -I $BOINCSRC/lib -DANDROID -c $SRC/src/main_nostos_boinc.cpp -o $W/main.o
$CXX $W/*.o $BOINC/api/libboinc_api.a $BOINC/lib/libboinc.a -static-libstdc++ -o $W/nostos_v${V}_android_aarch64
$NDK/llvm-strip $W/nostos_v${V}_android_aarch64
cp $W/nostos_v${V}_android_aarch64 $OUT/
$NDK/llvm-readelf -d $OUT/nostos_v${V}_android_aarch64 | grep NEEDED
ls -la $OUT/
