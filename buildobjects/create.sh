#!/bin/bash

# Add directory containing cl.exe to PATH
export PATH="$PATH:/Program Files/Microsoft Visual Studio/2022/Community/VC/Tools/MSVC/14.40.33807/bin/Hostx64/x64"

# Start compilation
echo "Starting compilation..."

cl.exe /Fo"accelerator.obj" /c /I../include /I"C:/Program Files/GSL/include" ../src/accelerator.cpp /D_USE_MATH_DEFINES /EHsc &
cl.exe /Fo"auxiliary.obj" /c /I../include /I"C:/Program Files/GSL/include" ../src/auxiliary.cpp /D_USE_MATH_DEFINES /EHsc &
cl.exe /Fo"commands.obj" /c /I../include /I"C:/Program Files/GSL/include" ../src/commands.cpp /D_USE_MATH_DEFINES /EHsc &
cl.exe /Fo"elements.obj" /c /I../include /I"C:/Program Files/GSL/include" ../src/elements.cpp /D_USE_MATH_DEFINES /EHsc &
cl.exe /Fo"exec.obj" /c /I../include /I"C:/Program Files/GSL/include" ../src/exec.cpp /D_USE_MATH_DEFINES /EHsc &
cl.exe /Fo"flat_file.obj" /c /I../include /I"C:/Program Files/GSL/include" ../src/flat_file.cpp /D_USE_MATH_DEFINES /EHsc &
cl.exe /Fo"kicktable.obj" /c /I../include /I"C:/Program Files/GSL/include" ../src/kicktable.cpp /D_USE_MATH_DEFINES /EHsc &
cl.exe /Fo"lattice.obj" /c /I../include /I"C:/Program Files/GSL/include" ../src/lattice.cpp /D_USE_MATH_DEFINES /EHsc &
cl.exe /Fo"linalg.obj" /c /I../include /I"C:/Program Files/GSL/include" ../src/linalg.cpp /D_USE_MATH_DEFINES /EHsc &
cl.exe /Fo"optics.obj" /c /I../include /I"C:/Program Files/GSL/include" ../src/optics.cpp /D_USE_MATH_DEFINES /EHsc &
cl.exe /Fo"output.obj" /c /I../include /I"C:/Program Files/GSL/include" ../src/output.cpp /D_USE_MATH_DEFINES /EHsc &
cl.exe /Fo"passmethods.obj" /c /I../include /I"C:/Program Files/GSL/include" ../src/passmethods.cpp /D_USE_MATH_DEFINES /EHsc &
cl.exe /Fo"tests.obj" /c /I../include /I"C:/Program Files/GSL/include" ../src/tests.cpp /D_USE_MATH_DEFINES /EHsc &
cl.exe /Fo"trackcpp.obj" /c /I../include /I"C:/Program Files/GSL/include" ../src/trackcpp.cpp /D_USE_MATH_DEFINES /EHsc &
cl.exe /Fo"tracking.obj" /c /I../include /I"C:/Program Files/GSL/include" ../src/tracking.cpp /D_USE_MATH_DEFINES /EHsc

lib.exe /OUT:libtrackcpp.lib \
    accelerator.obj \
    auxiliary.obj \
    commands.obj \
    elements.obj \
    exec.obj \
    flat_file.obj \
    kicktable.obj \
    lattice.obj \
    linalg.obj \
    optics.obj \
    output.obj \
    passmethods.obj \
    tests.obj \
    trackcpp.obj \
    tracking.obj \
    /link /LIBPATH:"C:\Program Files\GSL\lib" \
    gsl.lib gslcblas.lib


cl.exe \
    /Fe:trackcpp.exe \
    accelerator.obj \
    auxiliary.obj \
    commands.obj \
    elements.obj \
    exec.obj \
    flat_file.obj \
    kicktable.obj \
    lattice.obj \
    linalg.obj \
    optics.obj \
    output.obj \
    passmethods.obj \
    tests.obj \
    trackcpp.obj \
    tracking.obj \
    libtrackcpp.lib \
    /link /LIBPATH:"C:\Program Files\GSL\lib" \
    gsl.lib gslcblas.lib

swig.exe -c++ -python \
    -I..\\python_package -I..\\include \
    -I"C:\Program Files\GSL\include" \
    -D_USE_MATH_DEFINES \
    -o "pythonbuild\trackcpp_wrap.cxx" \
    "..\python_package\trackcpp.i"

cl.exe /Fo:"pythonbuild\interface.obj" /c ..\\python_package\\interface.cpp \
    /I..\\python_package /I..\\include \
    /I"C:\Program Files\GSL\include" \
    /D_USE_MATH_DEFINES

cl.exe /Fo:"pythonbuild\trackcpp_wrap.obj" /c pythonbuild\\trackcpp_wrap.cxx \
    /I..\\python_package /I..\\include \
    /I"C:\Program Files\GSL\include" \
    /I"C:\Users\vitin\miniforge3\envs\sirius\include" \
    /I"C:\Users\vitin\miniforge3\envs\sirius\Lib\site-packages\numpy\core\include" \
    /D_USE_MATH_DEFINES

cd pythonbuild


cl.exe /LD /Fo:"_tr.dll" "..\libtrackcpp.lib" trackcpp_wrap.obj interface.obj \
    /link /OUT:_tr.dll /LIBPATH:"C:\Program Files\GSL\lib" \
    /LIBPATH:"C:\Users\vitin\miniforge3\envs\sirius\libs" \
    gsl.lib gslcblas.lib python39.lib

# cl.exe /Fo:"_tr.obj" ..\\libtrackcpp.lib trackcpp_wrap.obj interface.obj \
#     /link /LIBPATH:"C:\Program Files\GSL\lib" \
#     /LIBPATH:"C:\Users\vitin\miniforge3\envs\sirius\libs" \
#     gsl.lib gslcblas.lib python39.lib

cd ..

echo "Finished!"

#cl /LD /Fo:"pythonbuild\_tr.pyd" /c libtrackcpp.lib interface.obj trackcpp_wrap.obj /link /LIBPATH:"C:\Program Files\GSL\lib" /LIBPATH:"C:\Users\vitin\miniforge3\envs\sirius\libs" gsl.lib gslcblas.lib python39.lib
