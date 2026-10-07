# Makefile for modulesLib Sep 01 13:09:26 2009

all: tideLib


PATH_OBJ = ../OBJ_INTEL
FFLAGS   = -w -O3 -xHost -qopt-zmm-usage=high -132 -traceback  -module $(PATH_OBJ)
FFLAGS7  = -w -O3 -xHost -qopt-zmm-usage=high -132 -traceback  -module $(PATH_OBJ)
FTN      = ifort
CFLAGS   = -I ./
CC       = icc


LIB = $(PATH_OBJ)/libtide.a

tideLib : $(LIB)

$(LIB) : \
	$(LIB)(ADMINT.o) \
        $(LIB)(ARG2.o) \
	$(LIB)(CAL2JD.o) \
	$(LIB)(DEHANTTIDEINEL.o) \
	$(LIB)(ETUTC.o) \
	$(LIB)(EVAL.o) \
	$(LIB)(HARDISP.o) \
	$(LIB)(JULDAT.o) \
	$(LIB)(LEAP.o) \
	$(LIB)(MDAY.o) \
	$(LIB)(NORM8.o) \
	$(LIB)(RECURS.o) \
	$(LIB)(SHELLS.o) \
	$(LIB)(SPLINE.o) \
	$(LIB)(SPROD.o) \
	$(LIB)(ST1IDIU.o) \
	$(LIB)(ST1ISEM.o) \
	$(LIB)(ST1L1.o) \
	$(LIB)(STEP2DIU.o) \
	$(LIB)(STEP2LON.o) \
	$(LIB)(TDFRPH2.o) \
	$(LIB)(TDFRPH.o) \
	$(LIB)(TOYMD.o) \
	$(LIB)(ZERO_VEC8.o) \
	$(LIB)(tide_com.o) \
	$(LIB)(tide_displace.o)

%.o : %.f90
	$(FTN) -c $(FFLAGS) $<
	ar -rv $(LIB) *.o
 
%.o : %.f
	$(FTN) -c $(FFLAGS7) $<
	ar -rv $(LIB) *.o

clean:
	rm -f $(LIB)
