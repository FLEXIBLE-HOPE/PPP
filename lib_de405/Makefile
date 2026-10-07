# Makefile for lib on SunOS 5.5.1 by unimake at Fri Dec  3 10:39:47 MET 1999

all: de405Lib


PATH_OBJ = ../OBJ_INTEL
FFLAGS   = -w -O3 -xHost -qopt-zmm-usage=high -132 -traceback  -module $(PATH_OBJ)
FFLAGS7  = -w -O3 -xHost -qopt-zmm-usage=high -132 -traceback  -module $(PATH_OBJ)
FTN      = ifort
CFLAGS   = -I ./
CC       = icc


LIB = $(PATH_OBJ)/libde405.a

de405Lib : $(LIB)

$(LIB) : \
        $(LIB)(interp.o) \
        $(LIB)(pleph.o) \
        $(LIB)(split.o) \
        $(LIB)(state.o)  

$(LIB)(state.o): jpleph.h


%.o : %.f90
	$(FTN) -c $(FFLAGS) $<
	ar -rv $(LIB) *.o
 
%.o : %.f
	$(FTN) -c $(FFLAGS7) $<
	ar -rv $(LIB) *.o

clean:
	rm -f $(LIB)
