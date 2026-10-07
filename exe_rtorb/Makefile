# Makefile for rtorb

.PHONY: all clean

all: rtorb


PATH_OBJ = ../OBJ_INTEL
FFLAGS   = -w -O3 -fp-model precise -132 -traceback -module $(PATH_OBJ) -mkl=parallel -qopt-zmm-usage=high
FTN      = ifort


DE405     = $(PATH_OBJ)/libde405.a
ITRS2GCRS = $(PATH_OBJ)/libitrs2gcrs.a
SOFA      = $(PATH_OBJ)/libsofa.a
COM       = $(PATH_OBJ)/libcom.a
MODULES   = $(PATH_OBJ)/libmodules.a

rtorb : $(DE405) $(ITRS2GCRS) $(SOFA) $(COM) $(MODULES)
	$(FTN) $(FFLAGS) rtorb.f90 get_rtorb_args.f90 $(DE405) $(ITRS2GCRS) $(SOFA) $(COM) $(MODULES) -o $(PATH_OBJ)/rtorb
	-rm -f *.o


clean:
	-rm -f $(PATH_OBJ)/rtorb
