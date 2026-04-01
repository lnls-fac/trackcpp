#### COMPILATION OPTIONS ####
PACKAGE = trackpybind
PREFIX ?=
PIP ?= pip
PYTHON ?= python3
ifeq ($(CONDA_PREFIX),)
	PREFIX = sudo -H
	PIP = pip-sirius
	PYTHON = python-sirius
endif

CXX = g++

CPPFLAGS = -O3 -Wall -shared -std=c++11 -fPIC
LIBS = $(shell gsl-config --libs)
LIBS += -ltrackcpp

INC = $(shell $(PYTHON) -m pybind11 --includes)
INC += $(shell gsl-config --cflags)
INC += -I../../include


.NOTPARALLEL:

all: $(PACKAGE)-install

$(PACKAGE)/trackcpp_pybind.so: 
	$(CXX) $(CPPFLAGS) $(INC) main.cpp $(LIBS) -o $@

$(PACKAGE)-pip: $(PACKAGE)/trackcpp_pybind.so
	$(PREFIX) $(PIP) install -e ./

$(PACKAGE)-stubs: $(PACKAGE)-pip
	$(shell pybind11-stubgen $(PACKAGE) && mv stubs/$(PACKAGE)/trackcpp_pybind.pyi stubs/$(PACKAGE)/__init__.pyi -t $(PACKAGE) && rm -rf stubs)

$(PACKAGE)-install: $(PACKAGE)-stubs

$(PACKAGE)-uninstall:
	$(PREFIX) $(PIP) uninstall -y $(PACKAGE)

clean:
	-rm -rf $(PACKAGE)/*.so $(PACKAGE)/*.pyi $(PACKAGE)/__pycache__ $(PACKAGE).egg-info
	
install: uninstall all

uninstall: clean $(PACKAGE)-uninstall

