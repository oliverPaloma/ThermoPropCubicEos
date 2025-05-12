ROOT_DIR:=$(dir $(abspath $(firstword $(MAKEFILE_LIST))))
C	= gcc
CC	= g++ #-Wall -pedantic
FC	= gfortran
#FCOPT	= -O2
FCOPT = -g -O0 
CCOPT	= -std=c++2a $(FCOPT)
ThermoProjPath	= /home/paloma/Documentos/ThermoPropCubicEos-main/ThermoPropCubicEos
ExecName = main
CPPFILES	= $(ExecName).cpp $(ThermoProjPath)/eos/eos.cpp \
			$(ThermoProjPath)/database/reading.cpp
HFILES	= $(ThermoProjPath)/includes.hpp $(ThermoProjPath)/eos/eos.hpp $(ThermoProjPath)/metodosMatematicos/metodosMatematicos.hpp \
			$(ThermoProjPath)/database/reading.hpp $(ThermoProjPath)/ThermoPropCubicEos.hpp
FFILES	= 
FORFILES	= 
CFILES		= 
OBJS	= $(CPPFILES:.cpp=.o) $(FFILES:.f=.o) $(FORFILES:.for=.o) $(CFILES:.c=.o)
EXE	= $(ExecName)

$(EXE)	: $(OBJS)
	  @$(CC) -o $(EXE) $(OBJS) $(GNUPLOT)

clean	:
	  rm *.o $(EXE)  $(ThermoProjPath)/eos/*.o $(ThermoProjPath)/database/*.o
	  
%.o:	%.cpp 
	$(CC) -c $(CCOPT) -o $@ $<

%.o:	%.cpp %.hpp
	$(CC) -c $(CCOPT) -o $@ $<

%.o:	%.for 
	$(FC) -c $(FCOPT) -o $@ $<

%.o:	%.f 
	$(FC) -c $(FCOPT) -o $@ $<

%.o:	%.c 
	$(C) -c $(FCOPT) -o $@ $<
