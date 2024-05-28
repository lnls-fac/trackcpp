#### COMPILATION OPTIONS ####
C = gcc
CXX = g++
SWIG = swig
MACHINE = -m64
OPT_FLAG = -O3
SOURCES_I = trackcpp.i
WRAPPEROBJS = trackcpp_wrap.o
INTERFACEOBJS = interface.o
TRACKCPPLIB = ..\build\libtrackcpp.a
PACKAGE = trackcpp
PREFIX = C:\Users\vitin\miniforge3\envs\sirius
PIP = $(PREFIX)\Scripts\pip.exe
PYTHON = $(PREFIX)\python.exe

PYTHON_INCLUDES := -I"$(PREFIX)\include"
NUMPY_INCLUDE := -I$(shell $(PYTHON) -c "import numpy; print(numpy.get_include())")
NUMPY_VERSION := $(shell $(PYTHON) -c "import numpy; print(numpy.__version__)")

OBJECTS = $(WRAPPEROBJS) $(INTERFACEOBJS)

CPPFLAGS = -std=c++11 -fPIC $(OPT_FLAG) -D_USE_MATH_DEFINES
LIBS = -L"$(PREFIX)\Library\lib" -lpthread -lgsl -lgslcblas -lm
LIBS += -L..\build -ltrackcpp
INC = -I"$(PREFIX)\Library\include"
INC += -I"$(PREFIX)\include"
INC += -I..\include $(PYTHON_INCLUDES) $(NUMPY_INCLUDE)

PYTHON_LIB_DIR := $(PREFIX)
PYTHON_LIB := python39

INTERFACESRCS = $(INTERFACEOBJS:.o=.cpp)
INTERFACEHDRS = $(INTERFACEOBJS:.o=.h)

.NOTPARALLEL:

all: trackcpp\_trackcpp.pyd

$(shell $(CXX) -MM $(INTERFACESRCS) > .depend)
-include .depend

trackcpp\_trackcpp.pyd: $(TRACKCPPLIB) $(OBJECTS)
	$(CXX) -shared -Wl,--whole-archive $(TRACKCPPLIB) -Wl,--no-whole-archive $(OBJECTS) $(LIBS) -L$(PYTHON_LIB_DIR) -l$(PYTHON_LIB) -o $@

trackcpp_wrap.cxx: $(SOURCES_I) $(INTERFACESRCS) $(INTERFACEHDRS)
	@if exist .\numpy.i ( \
	    echo "numpy.i already here, good" \
	) else ( \
	    curl -L -o numpy.i https://raw.githubusercontent.com/numpy/numpy/v$(NUMPY_VERSION)/tools/swig/numpy.i \
	)
	$(SWIG) -c++ -python $(INC) trackcpp.i
	@if exist trackcpp.py ( \
	    copy /Y trackcpp.py trackcpp \
	)

install: uninstall all
	$(PIP) install --no-deps ./

uninstall:
	$(PIP) uninstall -y $(PACKAGE)

develop-install: develop-uninstall all
	$(PIP) install --no-deps -e ./

# known issue: It will fail to uninstall scripts
#  if they were installed in develop mode
develop-uninstall:
	$(PIP) uninstall -y $(PACKAGE)

clean:
	@echo to be clean!
	if exist trackcpp\*.pyd del /Q /S trackcpp\*.pyd
	if exist *.o del /Q /S *.o
	if exist *_wrap.cxx del /Q /S *_wrap.cxx
	if exist trackcpp\__pycache__ del /Q /S trackcpp\__pycache__
	if exist trackcpp\trackcpp.py del /Q /S trackcpp\trackcpp.py
	if exist build rmdir /Q /S build
	if exist dist del /Q /S  dist
	if exist trackcpp.egg-info rmdir /Q /S trackcpp.egg-info

$(INTERFACESRCS): $(INTERFACEHDRS)

%.o: %.cpp
	$(CXX) -c $(CPPFLAGS) $(INC) $< -o $@

%.o: %.cxx
	$(CXX) -c $(CPPFLAGS) $(INC) $< -o $@
