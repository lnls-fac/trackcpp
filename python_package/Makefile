
#### COMPILATION OPTIONS ####
C = gcc
CXX = g++
MACHINE = -m64
OPT_FLAG = -O3
SOURCES_I =     trackcpp.i \
                elements.i
WRAPPEROBJS =   trackcpp_wrap.o
INTERFACEOBJS = interface.o \
                elementswrapper.o
TRACKCPPLIB = ../build/libtrackcpp.a

OBJECTS = $(WRAPPEROBJS) $(INTERFACEOBJS)

CPPFLAGS = -std=c++11 -fPIC $(OPT_FLAG)
LIBS = -lgsl -lblas -L../build -ltrackcpp
INC = -I/usr/include/python3.4 -I../include
PYTHON = python3
SETUPARGS =

INTERFACESRCS = $(INTERFACEOBJS:.o=.cpp)
INTERFACEHDRS = $(INTERFACEOBJS:.o=.h)

all: trackcpp/_trackcpp.so

$(shell $(CXX) -MM $(INTERFACESRCS) > .depend)
-include .depend

trackcpp/_trackcpp.so: $(TRACKCPPLIB) $(OBJECTS)
	$(CXX) -shared -Wl,--whole-archive $(TRACKCPPLIB) -Wl,--no-whole-archive $(OBJECTS) $(LIBS) -o trackcpp/_trackcpp.so

trackcpp_wrap.cxx: $(SOURCES_I) $(INTERFACESRCS) $(INTERFACEHDRS)
	swig -c++ -python $(INC) trackcpp.i && cp -f trackcpp.py trackcpp

install: all
	$(PYTHON) setup.py install $(SETUPARGS)

develop: all
	$(PYTHON) setup.py develop $(SETUPARGS)

clean:
	-rm -rf *.o trackcpp/*.so *_wrap.cxx trackcpp/__pycache__ trackcpp/trackcpp.py build dist trackcpp.egg-info

$(INTERFACESRCS): $(INTERFACEHDRS)


%.o: %.cpp
	$(CXX) -c $(CPPFLAGS) $(INC) $< -o $@

%.o: %.cxx
	$(CXX) -c $(CPPFLAGS) $(INC) $< -o $@
