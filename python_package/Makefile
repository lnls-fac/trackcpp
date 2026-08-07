
#### COMPILATION OPTIONS ####
C = gcc
CXX = g++
SWIG = swig
MACHINE = -m64
OPT_FLAG = -O3
SOURCES_I =     trackcpp.i
WRAPPEROBJS =   trackcpp_wrap.o
INTERFACEOBJS = interface.o
TRACKCPPLIB = ../build/libtrackcpp.a
PACKAGE = trackcpp
PREFIX ?=
PIP ?= pip
PYTHON ?= python3
ifeq ($(CONDA_PREFIX),)
	PREFIX = sudo -H
	PIP = pip-sirius
	PYTHON = python-sirius
endif

PYTHON_INCLUDES := $(shell $(PYTHON)-config --includes)
NUMPY_INCLUDE := -I$(shell $(PYTHON) -c "import numpy; print(numpy.get_include())")
NUMPY_VERSION := $(shell $(PYTHON) -c "import numpy; print(numpy.__version__)")

OBJECTS = $(WRAPPEROBJS) $(INTERFACEOBJS)

CPPFLAGS = -std=c++11 -fPIC $(OPT_FLAG)
LIBS = $(shell gsl-config --libs)
LIBS += -L../build -ltrackcpp
INC = $(shell gsl-config --cflags)
INC += -I../include $(PYTHON_INCLUDES) $(NUMPY_INCLUDE)

INTERFACESRCS = $(INTERFACEOBJS:.o=.cpp)
INTERFACEHDRS = $(INTERFACEOBJS:.o=.h)

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    SOFLAGS = -bundle -undefined dynamic_lookup
    ifneq ($(CONDA_PREFIX),)
        SOFLAGS += -Wl,-rpath,$(CONDA_PREFIX)/lib
    endif
    ARCHIVE_LINK = -Wl,-force_load,$(TRACKCPPLIB)
else
    SOFLAGS = -shared
    ARCHIVE_LINK = -Wl,--whole-archive $(TRACKCPPLIB) -Wl,--no-whole-archive
endif

.NOTPARALLEL:

all: trackcpp/_trackcpp.so

$(shell $(CXX) -MM $(INC) $(INTERFACESRCS) > .depend)
-include .depend

trackcpp/_trackcpp.so: $(TRACKCPPLIB) $(OBJECTS)
	$(CXX) $(SOFLAGS) $(ARCHIVE_LINK) $(OBJECTS) $(LIBS) -o trackcpp/_trackcpp.so

trackcpp_wrap.cxx: $(SOURCES_I) $(INTERFACESRCS) $(INTERFACEHDRS)
	[ -f ./numpy.i ] && echo "numpy.i already here, good" || \
	(command -v wget >/dev/null 2>&1 && wget https://raw.githubusercontent.com/numpy/numpy/v$(NUMPY_VERSION)/tools/swig/numpy.i || curl -sL -o numpy.i https://raw.githubusercontent.com/numpy/numpy/v$(NUMPY_VERSION)/tools/swig/numpy.i)
	$(SWIG) -c++ -python $(INC) trackcpp.i && cp -f trackcpp.py trackcpp

install: uninstall all
	$(PREFIX) $(PIP) install --no-deps ./

uninstall:
	$(PREFIX) $(PIP) uninstall -y $(PACKAGE)

develop-install: develop-uninstall all
	$(PIP) install --no-deps -e ./

# known issue: It will fail to uninstall scripts
#  if they were installed in develop mode
develop-uninstall:
	$(PIP) uninstall -y $(PACKAGE)

clean:
	-rm -rf *.o trackcpp/*.so *_wrap.cxx trackcpp/__pycache__ trackcpp/trackcpp.py build dist trackcpp.egg-info

$(INTERFACESRCS): $(INTERFACEHDRS)

%.o: %.cpp
	$(CXX) -c $(CPPFLAGS) $(INC) $< -o $@

%.o: %.cxx
	$(CXX) -c $(CPPFLAGS) $(INC) $< -o $@
