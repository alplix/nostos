@echo off
rem Nostos CUDA GPU app for Windows x86_64 (MSVC + nvcc, the same toolchain and BOINC MSVC libraries as
rem bitboinc's Keraunos CUDA build). Output: build\gpuwin\nostos_cuda_v%V%_windows_x86_64.exe
call "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat" >nul 2>nul
set NVCC="C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v13.3\bin\nvcc.exe"
set V=1
if not "%NOSTOS_VERSION%"=="" set V=%NOSTOS_VERSION%
set SRC=%~dp0..\src
set BOINC=C:\Users\Alp\Documents\bitboinc-app\gpu_test_local\boinc
set GENCODES=-gencode arch=compute_75,code=sm_75 -gencode arch=compute_80,code=sm_80 -gencode arch=compute_86,code=sm_86 -gencode arch=compute_89,code=sm_89 -gencode arch=compute_90,code=sm_90 -gencode arch=compute_100,code=sm_100 -gencode arch=compute_120,code=sm_120 -gencode arch=compute_120,code=compute_120
set OUT=%~dp0..\build\gpuwin
if not exist %OUT% mkdir %OUT%
cd /d %OUT%
for %%f in (coords orbital_elements potential montecarlo nostos_io) do cl /nologo /O2 /std:c11 /MT /I %SRC% /c %SRC%\%%f.c /Fo%%f.obj || exit /b 1
%NVCC% -O3 -std=c++17 --fmad=false %GENCODES% -I %SRC% -c %SRC%\nostos_cuda.cu -o nostos_cuda.obj -Xcompiler /MT || exit /b 1
%NVCC% -O2 -std=c++17 -DNG_CUDA -I %SRC% -I %BOINC%\api -I %BOINC%\lib -c %SRC%\main_nostos_gpu_boinc.cpp -o main.obj -Xcompiler /MT || exit /b 1
%NVCC% -o nostos_cuda_v%V%_windows_x86_64.exe main.obj nostos_cuda.obj coords.obj orbital_elements.obj potential.obj montecarlo.obj nostos_io.obj %BOINC%\win_build\Build\x64\Release\libboincapi.lib %BOINC%\win_build\Build\x64\Release\libboinc.lib -Xcompiler /MT -Xlinker advapi32.lib -Xlinker shell32.lib -Xlinker user32.lib -Xlinker ole32.lib -Xlinker psapi.lib -Xlinker /SUBSYSTEM:CONSOLE || exit /b 1
echo === DONE ===
