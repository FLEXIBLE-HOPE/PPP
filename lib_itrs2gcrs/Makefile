# Makefile for modulesLib Sep 01 13:09:26 2009

all: itrs2gcrsLib


PATH_OBJ = ../OBJ_INTEL
FFLAGS   = -w -O3 -xHost -qopt-zmm-usage=high -132 -traceback  -module $(PATH_OBJ)
FFLAGS7  = -w -O3 -xHost -qopt-zmm-usage=high -132 -traceback  -module $(PATH_OBJ)
FTN      = ifort
CFLAGS   = -I ./
CC       = icc


LIB = $(PATH_OBJ)/libitrs2gcrs.a

itrs2gcrsLib : $(LIB)

$(LIB) : \
	$(LIB)(CNMTX.o) \
        $(LIB)(FUNDARG.o) \
	$(LIB)(ORTHO_EOP.o) \
	$(LIB)(PMSDNUT2.o) \
	$(LIB)(RG_ZONT2.o) \
	$(LIB)(UTLIBR.o) \
	$(LIB)(itrs2gcrs.o) \
	$(LIB)(eop_partial.o) \
	$(LIB)(hfeop_xyu_mod.o) \
	$(LIB)(hf_eop.o)


%.o : %.f90
	$(FTN) -c $(FFLAGS) $<
	ar -rv $(LIB) *.o
 
%.o : %.f
	$(FTN) -c $(FFLAGS7) $<
	ar -rv $(LIB) *.o

clean:
	rm -f $(LIB)
