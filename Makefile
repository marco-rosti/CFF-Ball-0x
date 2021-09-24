
# gnu compiler
FC = gfortran
# FFLAGS := -O3 -ffixed-line-length-none -cpp -floop-nest-optimize #-mavx #-fno-math-errno
# FFLAGS := -O4 -fdefault-real-8 -ffixed-line-length-none -cpp -floop-nest-optimize -mtune=native -march=native #-fno-math-errno
FFLAGS := -O3 -fdefault-real-8 -ffixed-line-length-none -cpp  #-fno-math-errno
BIG  := -mcmodel=medium
DBG  := #-g -Wall -Wextra -fcheck=all -fbacktrace -fdump-core
PROF := #-pg
OMP  := -fopenmp
LIB  := #

# aocc (amd) compiler
# FC = flang
# FFLAGS := -O2 -Mextend -Mpreprocess #-fno-math-errno
# BIG  := -mcmodel=medium
# DBG  := #
# PROF := #
# OMP  := -fopenmp
# LIB  := #

# intel
# FC = ifort
# FFLAGS := -fpconstant -O3 -xHost -132 -cpp -fp-model consistent#-ipo
# BIG  := -mcmodel=medium -shared-intel
# DBG  := #-g -Wall -Wextra -fcheck=all#-traceback
# PROF := #-pg
# OMP  := -qopenmp
# LIB  := #

TARGET = ball0x.exe

SRC = param.f90 common.f90 stats.f90 init.f90 forcing.f90 timeIntegration.f90 destroy.f90 main.f90

OBJ = $(SRC:.f90=.o)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(FC) $(FFLAGS) $(BIG) $(DBG) $(PROF) $(OMP) -o $@ $(OBJ) $(LIB)

.PHONY: clean  
clean:
	rm -f *.o *.mod $(TARGET)

veryclean: 
	rm -f *.o *.mod $(TARGET) 

%.o: %.f90
	$(FC) $(INCLUDE) -c -o $@ $(FFLAGS) $(BIG) $(DBG) $(PROF) $(OMP) $<
common.o:          config.h param.f90
destroy.o:         config.h timeIntegration.f90 common.f90
forcing.o:         config.h common.f90 param.f90
identifyNNP.o:     forcing.f90 common.f90 param.f90
init.o:            config.h common.f90 param.f90
main.o:            destroy.f90 timeIntegration.f90 init.f90 param.f90
stats.o:           common.f90
timeIntegration.o: config.h forcing.o stats.o common.f90 param.f90
