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
VERSION=$(shell type ${FILE})

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
					linalg.cpp \
					auxiliary.cpp
BINSOURCES_CPP =	exec.cpp \
					tests.cpp \
					commands.cpp \


AUXFILES  = VERSION

LIBS = -L"C:\Users\vitin\miniforge3\envs\sirius\Library\lib"
LIBS += -lpthread -lgsl -lgslcblas -lm
INC = -I.\include
INC += -I"C:\Users\vitin\miniforge3\envs\sirius\Library\include"

PREFIX = "C:\Users\vitin\miniforge3\envs\sirius"

BINDEST_DIR = $(PREFIX)\Library\bin
LIBDEST_DIR = $(PREFIX)\Library\lib
INCDEST_DIR = $(PREFIX)\Library\include
OBJDIR = .\build
SRCDIR = .\src
INCDIR = .\include

PYTHON_PACKAGE_DIR = python_package

#### I STILL DONT GOT THIS
# $(shell type nul > $(SRCDIR)\output.cpp) # this is so that last compilation time always goes into executable

# WARNINGS_CFLAGS += -Wall
# WARNINGS_CFLAGS += -Wextra # reasonable and standard
# WARNINGS_CFLAGS += -Wshadow # warn the user if a variable declaration shadows one from a parent context
# WARNINGS_CFLAGS += -Wnon-virtual-dtor # warn the user if a class with virtual functions has a non-virtual destructor.
# #                                     # This helps catch hard to track down memory errors
# # WARNINGS_CFLAGS += -Wdouble-promotion # warn if float is implicit promoted to double
# # WARNINGS_CFLAGS += -Wformat=2 # warn on security issues around functions that format output (ie printf)
# # WARNINGS_CFLAGS += -Wold-style-cast # warn for c-style casts
# # WARNINGS_CFLAGS += -Wsign-conversion # warn on sign conversions
# WARNINGS_CFLAGS += -Wcast-align # warn for potential performance problem casts
# WARNINGS_CFLAGS += -Wconversion # warn on type conversions that may lose data
# # WARNINGS_CFLAGS += -Wduplicated-branches # warn if if \ else branches have duplicated code \\ not available in g++ 6.3.0
# WARNINGS_CFLAGS += -Wduplicated-cond # warn if if \ else chain has duplicated conditions
# # WARNINGS_CFLAGS += -Wimplicit-fallthrough # warn on statements that fallthrough without an explicit annotation  \\ not available in g++ 6.3.0
# WARNINGS_CFLAGS += -Wlogical-op # warn about logical operations being used where bitwise were probably wanted
# WARNINGS_CFLAGS += -Wmisleading-indentation # warn if indentation implies blocks where blocks do not exist
# WARNINGS_CFLAGS += -Wnull-dereference # warn if a null dereference is detected
# WARNINGS_CFLAGS += -Woverloaded-virtual # warn if you overload (not override) a virtual function
# WARNINGS_CFLAGS += -Wpedantic # warn if non-standard C++ is used
# WARNINGS_CFLAGS += -Wunused # warn on anything being unused
# WARNINGS_CFLAGS += -Wuseless-cast # warn if you perform a cast to the same type

ifeq ($(MAKECMDGOALS),trackcpp-debug)
  CFLAGS    = $(MACHINE) $(DBG_FLAG) $(DFLAGS) -pthread
else
  CFLAGS    = $(MACHINE) $(OPT_FLAG) $(DFLAGS) -pthread
endif
CFLAGS += $(WARNINGS_CFLAGS)
CFLAGS += -D_USE_MATH_DEFINES

LIBOBJECTS = $(addprefix $(OBJDIR)\, $(LIBSOURCES_CPP:.cpp=.o))
BINOBJECTS = $(addprefix $(OBJDIR)\, $(BINSOURCES_CPP:.cpp=.o))
LDFLAGS    = $(MACHINE)

.PHONY: all alllibs trackcpp clean cleanall

#### TARGETS ####

all:  libtrackcpp trackcpp python_package

# Rule to generate dependencies
depend: .depend

.depend:
	$(shell echo $(CXX) -MM $(CFLAGS) $(addprefix $(SRCDIR)\, $(LIBSOURCES_CPP)) $(addprefix $(SRCDIR)\, $(BINSOURCES_CPP)) > .depend.temp)
	$(shell type nul > .depend)
	@powershell -Command "(Get-Content .depend.temp) | %{ $_ -replace '\.o', '$(OBJDIR)\$$&' } | Set-Content .depend"
# @del .depend.temp

libtrackcpp: $(OBJDIR)\libtrackcpp.a

trackcpp: $(OBJDIR)\trackcpp

python_package: $(PYTHON_PACKAGE_DIR)\trackcpp\_trackcpp.so

$(PYTHON_PACKAGE_DIR)\trackcpp\_trackcpp.so: libtrackcpp
	$(MAKE) -C $(PYTHON_PACKAGE_DIR)

$(OBJDIR)\libtrackcpp.a: $(LIBOBJECTS)
	$(AR) $(ARFLAGS) $@ $^

$(OBJDIR)\trackcpp: libtrackcpp $(BINOBJECTS)
	$(CXX) $(LDFLAGS) $(BINOBJECTS) $(OBJDIR)\libtrackcpp.a $(LIBS) -o $@

$(LIBOBJECTS): | $(OBJDIR)

$(BINOBJECTS): | $(OBJDIR)

$(OBJDIR):
	mkdir $(OBJDIR)

install-cpp: uninstall-cpp libtrackcpp trackcpp
	if exist $(OBJDIR)\trackcpp.exe copy $(OBJDIR)\trackcpp.exe $(BINDEST_DIR)
	if exist $(OBJDIR)\libtrackcpp.a copy $(OBJDIR)\libtrackcpp.a $(LIBDEST_DIR)
#	if exist $(INCDIR)\trackcpp xcopy /E /I $(INCDIR)\trackcpp $(INCDEST_DIR)

uninstall-cpp:
	if exist $(BINDEST_DIR)\trackcpp.exe del /Q $(BINDEST_DIR)\trackcpp.exe
	if exist $(LIBDEST_DIR)\libtrackcpp.a del /Q $(LIBDEST_DIR)\libtrackcpp.a
# if exist $(INCDEST_DIR)\trackcpp del /Q  $(INCDEST_DIR)\trackcpp

install-py: uninstall-py python_package
	$(MAKE) install -C $(PYTHON_PACKAGE_DIR)

uninstall-py:
	$(MAKE) uninstall -C $(PYTHON_PACKAGE_DIR)

install: clean install-cpp install-py

uninstall: uninstall-cpp uninstall-py

clean:
	if exist $(OBJDIR) rmdir /S /Q $(OBJDIR)
	if exist trackcpp del /Q trackcpp
	if exist trackcpp-debug del /Q trackcpp-debug
# if exist .depend del /Q .depend
	if exist *.exe del /Q *.exe
	if exist *.dat del /Q *.dat
	if exist *~ del /Q *~
	if exist *.o del /Q *.o
	if exist *.a del /Q *.a

#### RULES ####

*.cpp: VERSION
	type nul > $(SRCDIR)\*.cpp
*.cc: VERSION
	type nul > $(SRCDIR)\*.cc
*.c: VERSION
	type nul > $(SRCDIR)\*.c

$(OBJDIR)\\%.o: $(SRCDIR)\%.cpp
	$(CXX) -c $(CFLAGS) $(INC) -I$(SRCDIR) -lm $< -o $@
