# General instructions

These general instructions will build the binary and library versions of
trackcpp, and also the python package. Before starting, read the requirements
in the sections below. If you intend

REQUIREMENTS

(See below).

INSTRUCTIONS

1. Compile with 'make all'.

1. Install with 'make install' or, alternatively, install a symbolic link to
   the local compiled versions with 'make develop'.

# trackcpp

REQUIREMENTS

- a C++ compiler with C++20 support (tested with GCC 4.8)
- pthread
- blas
- gsl (GNU Scientific Library)

For Debian systems the requirements can be installed using:

```command
# Debian
apt install libgsl-dev  libblas-dev swig
```

Clang tooling: As we are using Makefiles, we can use Bear to generate a compilation database file.

```command
apt install bear clang-tools clang-format clang-tidy
```

INSTRUCTIONS

1. Compile with 'make all'.

1. Install trackcpp with 'make install'. The installation directory is
   $(DEST_DIR), with DEST_DIR=/usr/local/bin by default. Alternatively, install a
   symbolic link to the local compiled version with 'make develop'.

# trackcpp Python package

REQUIREMENTS

- a C++ compiler with C++20 support (tested with GCC 4.8)
- blas
- gsl (GNU Scientific Library)
- SWIG 2 (tested with SWIG 2.0.11)
- Python 3 (tested with Python 3.4) with header files
- Python setuptools (>=16.0)
- libtrackcpp.a (see previous section)

INSTRUCTIONS

1. Check include and library paths in Makefile. On a Linux system, the default
   options should work (tested on Ubuntu 14.04 LTS); if necessary, alter the
   paths. An auxiliary script named 'compile' is provided and may be helpful on
   Mac OS X.

1. Run 'make all'.

1. Install the package with 'make install'. This will run 'python3 setup.py
   install', installing the Python package to a default location. Arguments to
   the command can be supplied through the variable SETUPARGS. Alternatively, the
   development version can be installed by running 'make develop', which will run
   'python3 setup.py develop' instead; arguments can be supplied the same way.
