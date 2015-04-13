
#### COMPILATION OPTIONS ####
C = gcc
CXX = g++
MACHINE = -m64
OPT_FLAG = -O3
SOURCES_I = trackcpp.i
INTERFACEOBJS = trackcpp_wrap.o \
		        interface.o \
		        elementswrapper.o
TRACKCPPOBJS = elements.o \
			   kicktable.o \
			   trackc++.o \
			   tracking.o \
			   passmethods.o \
			   lattice.o \
			   flat_file.o \
               optics.o

OBJECTS = $(TRACKCPPOBJS) $(INTERFACEOBJS)
CPPFLAGS = -std=c++11 -fPIC $(OPT_FLAG)
LIBS = -lgsl -lblas
INC = -I/usr/include/python3.4 -I../
prefix = ../../..

TRACKCPPSRCS = $(addprefix ../, $(TRACKCPPOBJS:.o=.cpp))

all: _trackcpp.so

$(shell $(CXX) -MM $(TRACKCPPSRCS) > .depend)
-include .depend

_trackcpp.so: $(OBJECTS)
	$(CXX) -shared $(OBJECTS) $(LIBS) -o _trackcpp.so

trackcpp_wrap.cxx: trackcpp.i
	swig -c++ -python trackcpp.i

install: all | $(prefix)/lib/python/trackcpp
	cp _trackcpp.so trackcpp.py __init__.py $(prefix)/lib/python/trackcpp

$(prefix)/lib/python/trackcpp:
	mkdir -p $(prefix)/lib/python/trackcpp

clean:
	-rm -rf *.o *.so *_wrap.cxx __pycache__ trackcpp.py


%.o: ../%.cpp
	$(CXX) -c $(CPPFLAGS) $(INC) $< -o $@

%.o: %.cpp
	$(CXX) -c $(CPPFLAGS) $(INC) $< -o $@

%.o: %.cxx
	$(CXX) -c $(CPPFLAGS) $(INC) $< -o $@
