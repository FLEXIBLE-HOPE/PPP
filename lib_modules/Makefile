# Makefile for modulesLib Sep 01 13:09:26 2009

all: modulesLib


PATH_OBJ = ../OBJ_INTEL
FFLAGS   = -w -O3 -132 -traceback  -module $(PATH_OBJ)
FFLAGS7  = -w -O3 -132 -traceback  -module $(PATH_OBJ)
FTN      = ifort
CFLAGS   = -I ./
CC       = icc


LIB = $(PATH_OBJ)/libmodules.a

modulesLib : $(LIB)

$(LIB) : \
	$(LIB)(mod_par.o) \
        $(LIB)(mod_const.o) \
	$(LIB)(mod_atx.o) \
	$(LIB)(mod_brdeph.o) \
	$(LIB)(mod_ckdctrl.o) \
	$(LIB)(mod_info.o) \
	$(LIB)(mod_ion.o) \
	$(LIB)(mod_obs.o) \
	$(LIB)(mod_orbit.o) \
	$(LIB)(mod_par.o) \
	$(LIB)(mod_sat.o) \
	$(LIB)(mod_sit.o) \
	$(LIB)(mod_crd.o) \
	$(LIB)(mod_tables.o) \
	$(LIB)(mod_amb.o) \
	$(LIB)(mod_qcctrl.o) \
	$(LIB)(mod_qcflag.o) 


%.o : %.f90
	$(FTN) -c $(FFLAGS) $<
	ar -rv $(LIB) *.o
 
%.o : %.f
	$(FTN) -c $(FFLAGS7) $<
	ar -rv $(LIB) *.o

clean:
	rm -f $(LIB) $(PATH_OBJ)/*.mod
