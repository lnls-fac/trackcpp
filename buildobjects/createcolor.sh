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
    cl.exe /O2 /Fo"objects\\${obj_file}" /c /I../include /I"C:/Program Files/GSL/include" "../src/${file}.cpp" /D_USE_MATH_DEFINES /EHsc
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
lib.exe /O2 /OUT:"objects\\libtrackcpp.lib" $(printf "objects\\%s.obj " "${sources[@]}") /LIBPATH:"C:\Program Files\GSL\lib" gsl.lib gslcblas.lib

wait

# Create executable
print_message $GREEN "Creating executable trackcpp.exe..."
cl.exe /O2 /Fe:"objects\\trackcpp.exe" $(printf "objects\\%s.obj " "${sources[@]}") "objects\\libtrackcpp.lib" /link /LIBPATH:"C:\Program Files\GSL\lib" gsl.lib gslcblas.lib /NODEFAULTLIB:LIBCMTD

wait

# Run SWIG
print_message $GREEN "Running SWIG..."
swig.exe -c++ -python -I..\\python_package -I..\\include -I"C:\Program Files\GSL\include" -D_USE_MATH_DEFINES -o "pythonbuild\trackcpp_wrap.cxx" "..\python_package\trackcpp.i"

wait

# Compile SWIG-generated wrapper and interface
print_message $CYAN "Compiling interface.cpp..."
cl.exe /O2 /Fo:"pythonbuild\interface.obj" /c ..\\python_package\\interface.cpp /I..\\python_package /I..\\include /I"C:\Program Files\GSL\include" /D_USE_MATH_DEFINES /EHsc

wait

print_message $CYAN "Compiling trackcpp_wrap.cxx..."
cl.exe /O2 /Fo:"pythonbuild\\trackcpp_wrap.obj" /c pythonbuild\\trackcpp_wrap.cxx /I..\\python_package /I..\\include /I"C:\Program Files\GSL\include" /I"C:\Users\vitin\miniforge3\envs\sirius\include" /I"C:\Users\vitin\miniforge3\envs\sirius\Lib\site-packages\numpy\core\include" /D_USE_MATH_DEFINES /EHsc

wait

cd pythonbuild
# Create Python extension
print_message $GREEN "Creating Python extension _tr.dll..."
cl.exe /O2 /LD /Fe:_trackcpp.pyd "..\\objects\\libtrackcpp.lib" trackcpp_wrap.obj interface.obj /link  /LIBPATH:"C:\Program Files\GSL\lib" /LIBPATH:"C:\Users\vitin\miniforge3\envs\sirius\libs" gsl.lib gslcblas.lib python39.lib /NODEFAULTLIB:LIBCMTD

wait

pwd
mv "_trackcpp.pyd" ../trackcpp
mv "trackcpp.py"  ../trackcpp

wait

cd ..
cmd.exe /c "pip.exe install -e ."

print_message $GREEN "Finished!"
