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

#### COMPILATION OPTIONS ####
CC          = gcc
CXX         = g++
AR          = ar
MACHINE     = -m64
OPT_FLAG    = -O3 -std=c++20 -fPIC
DBG_FLAG    = -O0 -g3 -std=c++20 -fPIC
ARFLAGS     = rcs
DFLAGS      = -DVERSION=$(VERSION)

AUXFILES  = VERSION

LIBS = -lgsl -lgslcblas -lpthread -lm
INC  = -I./include
BINDEST_DIR = /usr/local/bin
LIBDEST_DIR = /usr/local/lib
INCDEST_DIR = /usr/local/include

OBJDIR = build
SRCDIR = src
INCDIR = include

PYTHON_PACKAGE_DIR = python_package

ifeq ($(MAKECMDGOALS),trackcpp-debug)
  CFLAGS    = $(MACHINE) $(DBG_FLAG) $(DFLAGS) -pthread
else
  CFLAGS    = $(MACHINE) $(OPT_FLAG) $(DFLAGS) -pthread
endif

LDFLAGS    = $(MACHINE)


#### TARGETS ####

test:
	$(CXX) $(CFLAGS) $(INC) -I./$(SRCDIR) -c test-kickmap.cpp
	$(CXX) $(LDFLAGS) test-kickmap.o $(LIBDEST_DIR)/libtrackcpp.a $(LIBS) -o test-kickmap

#### RULES ####

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) -c $(CFLAGS) $(INC) -I./$(SRCDIR) $< -o $@;

clean:
	rm -rf *.o test-kickmap
