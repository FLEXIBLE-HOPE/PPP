# Makefile for lib on SunOS 5.5.1 by unimake at Fri Dec  3 10:39:47 MET 1999

all: lambdaLib

# ASSIGMENTS

PATH_OBJ = ../OBJ_INTEL
FFLAGS   = -w -O3 -xHost -qopt-zmm-usage=high -132 -traceback  -module $(PATH_OBJ)
FFLAGS7  = -w -O3 -xHost -qopt-zmm-usage=high -132 -traceback  -module $(PATH_OBJ)
FTN      = ifort
CFLAGS   = -I ./
CC       = icc

LIB = $(PATH_OBJ)/liblambda.a

lambdaLib : $(LIB)

$(LIB) : \
        $(LIB)(backts.o) \
        $(LIB)(chistrt4.o) \
        $(LIB)(collects.o) \
        $(LIB)(dinki.o) \
        $(LIB)(fi71.o) \
        $(LIB)(fmfac6.o) \
        $(LIB)(invlt2d.o) \
        $(LIB)(jnt2.o) \
        $(LIB)(lambda4.o) \
        $(LIB)(src1i.o) \
        $(LIB)(stores.o) \
        $(LIB)(ztrani.o) 
  
%.o : %.f90
	$(FTN) -c $(FFLAGS) $<
	ar -rv $(LIB) *.o
 
%.o : %.f
	$(FTN) -c $(FFLAGS7) $<
	ar -rv $(LIB) *.o

clean:
	rm -f $(LIB)
