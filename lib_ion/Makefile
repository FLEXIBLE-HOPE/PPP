# Makefile for modulesLib Sep 01 13:09:26 2009

all: ionLib


PATH_OBJ = ../OBJ_INTEL
FFLAGS   = -w -O3 -xHost -qopt-zmm-usage=high -132 -traceback  -module $(PATH_OBJ)
FFLAGS7  = -w -O3 -xHost -qopt-zmm-usage=high -132 -traceback  -module $(PATH_OBJ)
FTN      = ifort
CFLAGS   = -I ./
CC       = icc


LIB = $(PATH_OBJ)/libion.a
DE405     = $(PATH_OBJ)/libde405.a
ITRS2GCRS = $(PATH_OBJ)/libitrs2gcrs.a
SOFA      = $(PATH_OBJ)/libsofa.a
ION       = $(PATH_OBJ)/libion.a
TROPO     = $(PATH_OBJ)/libtropo.a
TIDE      = $(PATH_OBJ)/libtide.a
COM       = $(PATH_OBJ)/libcom.a
MODULES   = $(PATH_OBJ)/libmodules.a

ion : $(LIB) $(DE405) $(ITRS2GCRS) $(SOFA) $(TROPO) $(TIDE) $(COM) $(MODULES)
	$(FTN) $(FFLAGS) test_gim.f90 $(LIB) $(DE405) $(ITRS2GCRS) $(SOFA) $(TROPO) $(TIDE) $(COM) $(MODULES) -o ion

ionLib : $(LIB)

$(LIB) : \
	$(LIB)(klobuchar.o) \
	$(LIB)(ipp.o) \
	$(LIB)(gim.o) \
	$(LIB)(ion2.o) \
	$(LIB)(iono_corr_higher_order.o) \
	$(LIB)(igrf13.o) \
	$(LIB)(read_ionex.o) \

%.o : %.f90
	$(FTN) -c $(FFLAGS) $<
	ar -rv $(LIB) *.o
 
%.o : %.for
	$(FTN) -c $(FFLAGS7) $<
	ar -rv $(LIB) *.o

clean:
	rm -f $(LIB)
