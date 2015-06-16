# TRACKC++
# ========
# Author: 		Ximenes R. Resende
# email:  		xresende@gmail.com, ximenes.resende@lnls.br
# affiliation:	LNLS - Laboratorio Nacional de Luz Sincrotron
# Date: 		Tue Dec 10 17:57:20 BRST 2013

#### READS LIB VERSION ####

FILE=VERSION
VERSION=$(shell cat ${FILE})

$(shell touch output.cpp) # this is so that last compilation time always goes into executable

#### COMPILATION OPTIONS ####
CC		    = gcc
CXX		    = g++
MACHINE		= -m64
OPT_FLAG	= -O3 -std=c++11
DBG_FLAG	= -O0 -g3 -std=c++11
DFLAGS      = -DVERSION=$(VERSION)
SOURCES_C   =
SOURCES_CPP	= 	lattice.cpp \
			elements.cpp \
			passmethods.cpp \
			tracking.cpp \
			trackcpp.cpp \
			tests.cpp \
			sirius_v500.cpp \
			flat_file.cpp optics.cpp \
			dynap.cpp commands.cpp \
			output.cpp \
			kicktable.cpp \
			exec.cpp
AUXFILES  = VERSION

LIBS            = tracking_mp/build/tracking_mp.a -lpthread -lgsl -lgslcblas
INC             =
DEST_DIR = /usr/local/bin
OBJDIR = build

ifeq ($(MAKECMDGOALS),trackcpp-debug)
	CFLAGS		= $(MACHINE) $(DBG_FLAG) $(DFLAGS)
else
	CFLAGS		= $(MACHINE) $(OPT_FLAG) $(DFLAGS)
endif

OBJECTS		= $(addprefix $(OBJDIR)/, $(SOURCES_CPP:.cpp=.o) $(SOURCES_C:.c=.o))
LDFLAGS		= $(MACHINE)

#### DERIVED CONDITIONALS AND VARIABLES ####
ifeq ($(shell hostname), uv100)
    CC          = /progs/users/gcc-4.7.3/bin/gcc
    CXX         = /progs/users/gcc-4.7.3/bin/g++
    LIBS        = -Wl,-rpath,/progs/users/gcc-4.7.3/lib64 tracking_mp/build/tracking_mp.a -lpthread ~/bin/GSL/lib/libgsl.a ~/bin/GSL/lib/libgslcblas.a
    INC         = -I../bin/GSL/include/
endif


.PHONY: all alllibs trackcpp clean cleanall

#### TARGETS ####

all:	alllibs trackcpp

#### GENERATES DEPENDENCY FILE ####
$(shell $(CXX) -MM $(SOURCES_CPP) $(SOURCES_C) | sed 's/.*\.o/$(OBJDIR)\/&/' > .depend)
-include .depend

alllibs:
	cd tracking_mp; make all;

trackcpp: $(OBJDIR)/trackcpp

trackcpp-debug:	$(OBJECTS)
	$(CXX) $(LDFLAGS) $(OBJECTS) $(LIBS) -o $(OBJDIR)/$@

$(OBJDIR)/trackcpp: $(OBJECTS)
	$(CXX) $(LDFLAGS) $(OBJECTS) $(LIBS) -o $@

$(OBJECTS): | $(OBJDIR)

$(OBJDIR):
	mkdir $(OBJDIR)

install: uninstall all | $(DEST_DIR)
	cp $(OBJDIR)/trackcpp $(DEST_DIR)

develop: uninstall all
	ln -srf $(OBJDIR)/trackcpp $(DEST_DIR)

$(DEST_DIR):
	mkdir $(DEST_DIR)

clean:
	-rm -rf $(OBJDIR) trackcpp trackcpp-debug .depend *.out *.dat *~

uninstall:
	-rm -rf $(DEST_DIR)/trackcpp

cleanall: clean
	cd tracking_mp; make clean;


#### RULES ####

*.cpp: VERSION
	touch *.cpp
*.cc: VERSION
		touch *.cc
*.c: VERSION
		touch *.c

$(OBJDIR)/%.o: %.c
	$(CC) -c $(CFLAGS) $(INC) $< -o $@

$(OBJDIR)/%.o: %.cc
	$(CXX) -c $(CFLAGS) $(INC) $< -o $@

$(OBJDIR)/%.o: %.cpp
	$(CXX) -c $(CFLAGS) $(INC) $< -o $@
