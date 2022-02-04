## TRACKCPP
## ========
## Author:      Accelerator Physics Group - LNLS
## contact:     xresende@gmail.com
## affiliation: Laboratorio Nacional de Luz Sincrotron
##
## The MIT License (MIT)
##
## Copyright (c) <year> <copyright holders>
##
## Permission is hereby granted, free of charge, to any person obtaining a copy
## of this software and associated documentation files (the "Software"), to deal
## in the Software without restriction, including without limitation the rights
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
## copies of the Software, and to permit persons to whom the Software is
## furnished to do so, subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in
## all copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
## AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
## OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
## THE SOFTWARE.

#### READS LIB VERSION ####

FILE=VERSION
VERSION=$(shell cat ${FILE})

#### COMPILATION OPTIONS ####
CC          = gcc
CXX         = g++
AR          = ar
MACHINE     = -m64
OPT_FLAG    = -O3 -std=c++11 -fPIC
DBG_FLAG    = -O0 -g3 -std=c++11 -fPIC
ARFLAGS     = rcs
DFLAGS      = -DVERSION=$(VERSION)
LIBSOURCES_CPP  =	lattice.cpp \
					elements.cpp \
					passmethods.cpp \
					tracking.cpp \
					trackcpp.cpp \
					flat_file.cpp \
					optics.cpp \
					diffusion_matrix.cpp \
					dynap.cpp \
					output.cpp \
					kicktable.cpp \
					multithreads.cpp \
					accelerator.cpp \
					naff.cpp \
					linalg.cpp
BINSOURCES_CPP =	exec.cpp \
					tests.cpp \
					commands.cpp \


AUXFILES  = VERSION

LIBS = -lgsl -lgslcblas -lpthread -lm
INC  = -I./include

ifeq ($(CONDA_PREFIX),)
    PREFIX = /usr/local
else
    PREFIX = $(CONDA_PREFIX)
endif
BINDEST_DIR = $(PREFIX)/bin
LIBDEST_DIR = $(PREFIX)/lib
INCDEST_DIR = $(PREFIX)/include

OBJDIR = build
SRCDIR = src
INCDIR = include

PYTHON_PACKAGE_DIR = python_package

$(shell touch $(SRCDIR)/output.cpp) # this is so that last compilation time always goes into executable

ifeq ($(MAKECMDGOALS),trackcpp-debug)
  CFLAGS    = $(MACHINE) $(DBG_FLAG) $(DFLAGS) -pthread
else
  CFLAGS    = $(MACHINE) $(OPT_FLAG) $(DFLAGS) -pthread
endif

LIBOBJECTS  = $(addprefix $(OBJDIR)/, $(LIBSOURCES_CPP:.cpp=.o))
BINOBJECTS  = $(addprefix $(OBJDIR)/, $(BINSOURCES_CPP:.cpp=.o))
LDFLAGS    = $(MACHINE)

#### DERIVED CONDITIONALS AND VARIABLES ####
ifeq ($(shell hostname), uv100)
    CC          = /progs/users/gcc-4.7.3/bin/gcc
    CXX         = /progs/users/gcc-4.7.3/bin/g++
    LIBS        = -Wl,-rpath,/progs/users/gcc-4.7.3/lib64 ~/bin/GSL/lib/libgsl.a ~/bin/GSL/lib/libgslcblas.a -lpthread
    INC         = -I../bin/GSL/include/
endif


.PHONY: all alllibs trackcpp clean cleanall

#### TARGETS ####

all:  libtrackcpp trackcpp python_package

#### GENERATES DEPENDENCY FILE ####
$(shell $(CXX) -MM $(CFLAGS) $(addprefix $(SRCDIR)/, $(LIBSOURCES_CPP)) $(addprefix $(SRCDIR)/, $(BINSOURCES_CPP)) | sed 's/.*\.o/$(OBJDIR)\/&/' > .depend)
-include .depend

libtrackcpp: $(OBJDIR)/libtrackcpp.a

trackcpp: $(OBJDIR)/trackcpp

python_package: $(PYTHON_PACKAGE_DIR)/trackcpp/_trackcpp.so

$(PYTHON_PACKAGE_DIR)/trackcpp/_trackcpp.so: libtrackcpp
	$(MAKE) -C $(PYTHON_PACKAGE_DIR)

$(OBJDIR)/libtrackcpp.a: $(LIBOBJECTS)
	$(AR) $(ARFLAGS) $@ $^

$(OBJDIR)/trackcpp: libtrackcpp $(BINOBJECTS)
	$(CXX) $(LDFLAGS) $(BINOBJECTS) $(OBJDIR)/libtrackcpp.a $(LIBS) -o $@

$(LIBOBJECTS): | $(OBJDIR)

$(BINOBJECTS): | $(OBJDIR)

$(OBJDIR):
	mkdir $(OBJDIR)

install-cpp: uninstall-cpp all
	cp $(OBJDIR)/trackcpp $(BINDEST_DIR)
	cp $(OBJDIR)/libtrackcpp.a $(LIBDEST_DIR)
	cp -r $(INCDIR)/trackcpp $(INCDEST_DIR)

uninstall-cpp:
	-rm -rf $(BINDEST_DIR)/trackcpp
	-rm -rf $(LIBDEST_DIR)/libtrackcpp.a
	-rm -rf $(INCDEST_DIR)/trackcpp

install-py: uninstall-py
	$(MAKE) install -C $(PYTHON_PACKAGE_DIR)

uninstall-py:
	$(MAKE) uninstall -C $(PYTHON_PACKAGE_DIR)

install: install-cpp install-py

uninstall: uninstall-cpp uninstall-py

develop-install-py: develop-uninstall-py all
	$(MAKE) develop-install -C $(PYTHON_PACKAGE_DIR)

develop-uninstall-py:
	$(MAKE) develop-uninstall -C $(PYTHON_PACKAGE_DIR)

$(BINDEST_DIR):
	mkdir $(BINDEST_DIR)

$(LIBDEST_DIR):
	mkdir $(LIBDEST_DIR)

$(INCDEST_DIR):
	mkdir $(INCDEST_DIR)

clean:
	-rm -rf $(OBJDIR) trackcpp trackcpp-debug .depend *.out *.dat *~ *.o *.a
	$(MAKE) clean -C $(PYTHON_PACKAGE_DIR)

cleanall: clean
	cd tracking_mp; make clean;


#### RULES ####

*.cpp: VERSION
	touch $(SRCDIR)/*.cpp
*.cc: VERSION
	touch $(SRCDIR)/*.cc
*.c: VERSION
	touch $(SRCDIR)/*.c

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) -c $(CFLAGS) $(INC) -I./$(SRCDIR) $< -o $@;
