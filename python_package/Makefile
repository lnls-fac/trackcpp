
#### COMPILATION OPTIONS ####
C = gcc
CXX = g++
MACHINE = -m64
OPT_FLAG = -O3
SOURCES_I = trackcpp.i \
			elements.i
WRAPPEROBJS = trackcpp_wrap.o
INTERFACEOBJS = interface.o \
		        elementswrapper.o
TRACKCPPOBJS = elements.o \
			   kicktable.o \
			   trackcpp.o \
			   tracking.o \
			   passmethods.o \
			   lattice.o \
			   flat_file.o \
               optics.o

OBJECTS = $(WRAPPEROBJS) $(TRACKCPPOBJS) $(INTERFACEOBJS)
CPPFLAGS = -std=c++11 -fPIC $(OPT_FLAG)
LIBS = -lgsl -lblas
INC = -I/usr/include/python3.4 -I../
PYTHON = python3
SETUPARGS =

TRACKCPPSRCS = $(addprefix ../, $(TRACKCPPOBJS:.o=.cpp))
INTERFACESRCS = $(INTERFACEOBJS:.o=.cpp)
INTERFACEHDRS = $(INTERFACEOBJS:.o=.h)

all: trackcpp/_trackcpp.so

$(shell $(CXX) -MM $(TRACKCPPSRCS) $(INTERFACESRCS) > .depend)
-include .depend

trackcpp/_trackcpp.so: $(OBJECTS)
	$(CXX) -shared $(OBJECTS) $(LIBS) -o trackcpp/_trackcpp.so

trackcpp_wrap.cxx: $(SOURCES_I) $(TRACKCPPSRCS) $(INTERFACESRCS) $(INTERFACEHDRS)
	swig -c++ -python trackcpp.i && cp -f trackcpp.py trackcpp

install: all
	$(PYTHON) setup.py install $(SETUPARGS)

develop: all
	$(PYTHON) setup.py develop $(SETUPARGS)

clean:
	-rm -rf *.o trackcpp/*.so *_wrap.cxx trackcpp/__pycache__ trackcpp/trackcpp.py build dist trackcpp.egg-info

$(INTERFACESRCS): $(INTERFACEHDRS)


%.o: ../%.cpp
	$(CXX) -c $(CPPFLAGS) $(INC) $< -o $@

%.o: %.cpp
	$(CXX) -c $(CPPFLAGS) $(INC) $< -o $@

%.o: %.cxx
	$(CXX) -c $(CPPFLAGS) $(INC) $< -o $@