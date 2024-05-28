#!/bin/bash

# Function to print colored messages
print_message() {
    local color_code="$1"
    local message="$2"
    echo -e "\n\033[${color_code}m${message}\033[0m"
}

# Define color codes
RED="31"
GREEN="32"
YELLOW="33"
BLUE="34"
MAGENTA="35"
CYAN="36"

# Add directory containing cl.exe to PATH
# export PATH="$PATH;/Program Files/Microsoft Visual Studio/2022/Community/VC/Auxiliary/Build"
# cmd.exe /C vcvarsall.bat x64

bash clean.sh

export PATH="$PATH;/Program Files/Microsoft Visual Studio/2022/Community/VC/Tools/MSVC/14.40.33807/bin/Hostx64/x64"
# export PATH="$PATH:C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.40.33807\bin\Hostx64\x64"

# Start compilation
print_message $GREEN "Starting compilation..."

# Compile each source file
compile() {
    local file="$1"
    local obj_file="$2"
    print_message $CYAN "Compiling ${file}..."
    cl.exe /Fo"${obj_file}" /c /I../include /I"C:/Program Files/GSL/include" "../src/${file}.cpp" /D_USE_MATH_DEFINES /EHsc
}

# List of source files
sources=("accelerator" "auxiliary" "commands" "elements" "exec" "flat_file" "kicktable" "lattice" "linalg" "optics" "output" "passmethods" "tests" "trackcpp" "tracking")

# Compile all source files
for src in "${sources[@]}"; do
    compile $src "${src}.obj"
done

wait

# Create library
print_message $GREEN "Creating library libtrackcpp.lib..."
lib.exe /OUT:libtrackcpp.lib $(printf "%s.obj " "${sources[@]}") /LIBPATH:"C:\Program Files\GSL\lib" gsl.lib gslcblas.lib

# Create executable
print_message $GREEN "Creating executable trackcpp.exe..."
cl.exe /Fe:trackcpp.exe $(printf "%s.obj " "${sources[@]}") libtrackcpp.lib /link /LIBPATH:"C:\Program Files\GSL\lib" gsl.lib gslcblas.lib /NODEFAULTLIB:LIBCMTD

# Run SWIG
print_message $GREEN "Running SWIG..."
swig.exe -c++ -python -I..\\python_package -I..\\include -I"C:\Program Files\GSL\include" -D_USE_MATH_DEFINES -o "pythonbuild\trackcpp_wrap.cxx" "..\python_package\trackcpp.i"

# Compile SWIG-generated wrapper and interface
print_message $CYAN "Compiling interface.cpp..."
cl.exe /Fo:"pythonbuild\interface.obj" /c ..\\python_package\\interface.cpp /I..\\python_package /I..\\include /I"C:\Program Files\GSL\include" /D_USE_MATH_DEFINES

print_message $CYAN "Compiling trackcpp_wrap.cxx..."
cl.exe /Fo:"pythonbuild\trackcpp_wrap.obj" /c pythonbuild\\trackcpp_wrap.cxx /I..\\python_package /I..\\include /I"C:\Program Files\GSL\include" /I"C:\Users\vitin\miniforge3\envs\sirius\include" /I"C:\Users\vitin\miniforge3\envs\sirius\Lib\site-packages\numpy\core\include" /D_USE_MATH_DEFINES

cd pythonbuild

# Create Python extension
print_message $GREEN "Creating Python extension _tr.dll..."
# cl.exe /LD /Fo:"_tr.dll" "..\libtrackcpp.lib" trackcpp_wrap.obj interface.obj /link /OUT:_tr.dll /LIBPATH:"C:\Program Files\GSL\lib" /LIBPATH:"C:\Users\vitin\miniforge3\envs\sirius\libs" gsl.lib gslcblas.lib python39.lib
# cl.exe /LD "..\libtrackcpp.lib" trackcpp_wrap.obj interface.obj /link  /LIBPATH:"C:\Program Files\GSL\lib" /LIBPATH:"C:\Users\vitin\miniforge3\envs\sirius\libs" gsl.lib gslcblas.lib python39.lib /NODEFAULTLIB:LIBCMTD
# cl.exe /D_USRDLL /D_WINDLL ..\\libtrackcpp.lib trackcpp_wrap.obj interface.obj /MT /link /LIBPATH:"C:\Program Files\GSL\lib" /LIBPATH:"C:\Users\vitin\miniforge3\envs\sirius\libs" gsl.lib gslcblas.lib python39.lib /DLL /OUT:testDLL.dll /NODEFAULTLIB:LIBCMTD
# cl.exe /LD "..\libtrackcpp.lib" trackcpp_wrap.obj interface.obj /link  /LIBPATH:"C:\Program Files\GSL\lib" /LIBPATH:"C:\Users\vitin\miniforge3\envs\sirius\libs" gsl.lib gslcblas.lib python39.lib /NODEFAULTLIB:LIBCMTD
cl.exe /LD /Fe:_trackcpp.pyd "..\libtrackcpp.lib" trackcpp_wrap.obj interface.obj /link  /LIBPATH:"C:\Program Files\GSL\lib" /LIBPATH:"C:\Users\vitin\miniforge3\envs\sirius\libs" gsl.lib gslcblas.lib python39.lib /NODEFAULTLIB:LIBCMTD

cd ..

print_message $GREEN "Finished!"
