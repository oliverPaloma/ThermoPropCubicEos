ROOT_DIR:=$(dir $(abspath $(firstword $(MAKEFILE_LIST))))
AutodiffPath    = $(ROOT_DIR)../../external/autodiff
EigenExternal   = $(ROOT_DIR)../../external/eigen
EigenSystem     = /usr/include/eigen3

ifdef EIGEN_PATH
EigenPath       = $(EIGEN_PATH)
else ifneq (,$(wildcard $(EigenExternal)/Eigen/Core))
EigenPath       = $(EigenExternal)
else
EigenPath       = $(EigenSystem)
endif
C       = gcc
CC      = g++ #-Wall -pedantic
FC      = gfortran
#FCOPT  = -O2
FCOPT = -g -O0
CCOPT   = -std=c++2a $(FCOPT) -DAUTODIFF_EIGEN_FOUND -I$(AutodiffPath) -I$(EigenPath)
ThermoProjPath  = $(ROOT_DIR)../..
ExecName = main
CPPFILES        = $(ExecName).cpp $(ThermoProjPath)/eos/eos.cpp \
                        $(ThermoProjPath)/database/reading.cpp $(ThermoProjPath)/equilibrium/equilibrium.cpp
HFILES  = $(ThermoProjPath)/includes.hpp $(ThermoProjPath)/eos/eos.hpp $(ThermoProjPath)/metodosMatematicos/metodosMatematicos.hpp \
                        $(ThermoProjPath)/database/reading.hpp $(ThermoProjPath)/ThermoPropCubicEos.hpp
FFILES  =
FORFILES        =
CFILES          =
OBJS    = $(CPPFILES:.cpp=.o) $(FFILES:.f=.o) $(FORFILES:.for=.o) $(CFILES:.c=.o)
EXE     = $(ExecName)

$(EXE)  : $(OBJS)
	@$(CC) -o $(EXE) $(OBJS) -lgfortran

clean   :
	rm -f *.o $(EXE)  $(ThermoProjPath)/eos/*.o $(ThermoProjPath)/database/*.o $(ThermoProjPath)/equilibrium/*.o

%.o:    %.cpp
	$(CC) -c $(CCOPT) -o $@ $<

%.o:    %.cpp %.hpp
	$(CC) -c $(CCOPT) -o $@ $<

%.o:    %.for
	$(FC) -c $(FCOPT) -o $@ $<

%.o:    %.f
	$(FC) -c $(FCOPT) -o $@ $<

%.o:    %.c
	$(C) -c $(FCOPT) -o $@ $<
