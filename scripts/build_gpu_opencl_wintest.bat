@echo off
rem Local test build of the Nostos OpenCL app (MSVC, NVIDIA vendor filter) -- release builds are made
rem on the bitboinc server with mingw (scripts/build_gpu_server.sh).
call "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat" >nul 2>nul
set SRC=%~dp0..\src
set BOINC=C:\Users\Alp\Documents\bitboinc-app\gpu_test_local\boinc
set OUT=%~dp0..\build\clwin
if not exist %OUT% mkdir %OUT%
cd /d %OUT%
for %%f in (coords orbital_elements potential montecarlo nostos_io) do cl /nologo /O2 /std:c11 /MT /I %SRC% /c %SRC%\%%f.c /Fo%%f.obj || exit /b 1
cl /nologo /O2 /std:c11 /MT /I %SRC% /I C:\temp /DCL_TARGET_OPENCL_VERSION=120 /DEXPECTED_GPU_VENDOR=\"nvidia\" /c %SRC%\nostos_opencl.c /Fonostos_opencl.obj || exit /b 1
cl /nologo /O2 /std:c11 /MT /I %SRC% /I C:\temp /DCL_TARGET_OPENCL_VERSION=120 /c %SRC%\cl_dynload.c /Focl_dynload.obj || exit /b 1
cl /nologo /O2 /std:c++17 /EHsc /MT /DNG_OPENCL /I %SRC% /I %BOINC%\api /I %BOINC%\lib /c %SRC%\main_nostos_gpu_boinc.cpp /Fomain.obj || exit /b 1
link /nologo main.obj nostos_opencl.obj cl_dynload.obj coords.obj orbital_elements.obj potential.obj montecarlo.obj nostos_io.obj %BOINC%\win_build\Build\x64\Release\libboincapi.lib %BOINC%\win_build\Build\x64\Release\libboinc.lib advapi32.lib shell32.lib user32.lib ole32.lib psapi.lib /OUT:nostos_opencl_test_nvidia.exe || exit /b 1
echo === DONE ===
