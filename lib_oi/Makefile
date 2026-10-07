# Makefile for oi on SunOS 5.5.1 by unimake at Fri Nov 20 13:09:26 MET 1998

all: oi


PATH_OBJ = ../OBJ_INTEL
FFLAGS   = -w -O3 -xHost -qopt-zmm-usage=high -132 -traceback  -module $(PATH_OBJ)
FFLAGS7  = -w -O3 -xHost -qopt-zmm-usage=high -132 -traceback  -module $(PATH_OBJ)
FTN      = ifort
CFLAGS   = -I ./
CC       = icc


LIB       = liboi.a
DE405     = $(PATH_OBJ)/libde405.a
ITRS2GCRS = $(PATH_OBJ)/libitrs2gcrs.a
SOFA      = $(PATH_OBJ)/libsofa.a
COM       = $(PATH_OBJ)/libcom.a
MODULES   = $(PATH_OBJ)/libmodules.a

oi : $(LIB) $(DE405) $(ITRS2GCRS) $(SOFA) $(COM) $(MODULES)
	$(FTN) $(FFLAGS) oi.f90 -Wl,--start-group $(LIB) -Wl,--end-group $(DE405) $(ITRS2GCRS) $(SOFA) $(COM) $(MODULES) -o $(PATH_OBJ)/oi
	-rm -f *.o
	mv $(LIB) $(PATH_OBJ)

$(LIB) : \
	$(LIB)(cprn2sys.o) \
	$(LIB)(read_thrust.o) \
	$(LIB)(dtm94.o) \
	$(LIB)(dtm2012.o) \
	$(LIB)(doodson_arg.o) \
	$(LIB)(eci2orb.o) \
	$(LIB)(get_oi_args.o) \
	$(LIB)(normalize_coeff.o) \
	$(LIB)(oi_adams_coef.o) \
	$(LIB)(oi_adams.o) \
	$(LIB)(oi_aod.o) \
	$(LIB)(oi_atm_drag.o) \
	$(LIB)(oi_customer_acc.o) \
	$(LIB)(oi_eph_planet.o) \
	$(LIB)(oi_etpole_tide.o) \
	$(LIB)(oi_fright.o) \
	$(LIB)(oi_gravity_pines.o) \
	$(LIB)(oi_meg_flnorb.o) \
	$(LIB)(oi_fes2004_tide.o) \
	$(LIB)(oi_fes2014b_tide.o) \
	$(LIB)(oi_eot11a_tide.o) \
	$(LIB)(oi_ocpole_tide.o) \
	$(LIB)(oi_planet_init.o) \
	$(LIB)(oi_point_mass.o) \
	$(LIB)(oi_relativity.o) \
	$(LIB)(oi_rkf_coef.o) \
	$(LIB)(oi_rkf.o) \
	$(LIB)(oi_shadow_factor.o) \
	$(LIB)(oi_solid_tides.o) \
	$(LIB)(oi_srp_bern.o) \
	$(LIB)(oi_srp_boxw.o) \
	$(LIB)(oi_srp_pri.o) \
	$(LIB)(oi_srp_gls.o) \
	$(LIB)(oi_srp_bds.o) \
	$(LIB)(oi_srp_gen.o) \
	$(LIB)(oi_thrust_acc.o) \
	$(LIB)(oi_trp_boxw.o) \
	$(LIB)(oi_sat_panel.o) \
	$(LIB)(read_ics.o) \
	$(LIB)(read_tide_freq.o) \
	$(LIB)(scf2ecef.o) \
	$(LIB)(SRPFBOXW.o) \
	$(LIB)(ERPFBOXW.o) \
	$(LIB)(PROPBOXW.o) \
	$(LIB)(SURFBOXW.o) \
        $(LIB)(getAcc.o) \
        $(LIB)(myinterpolation.o) \
        $(LIB)(radiationForce.o) \
        $(LIB)(uclsrp.o)  \
        $(LIB)(ucl_grid_reader.o) \


%.o : %.f90
	$(FTN) -c $(FFLAGS) $<
	ar -rv $(LIB) *.o
 
%.o : %.f
	$(FTN) -c $(FFLAGS7) $<
	ar -rv $(LIB) *.o

clean:
	rm -f $(LIB) $(PATH_OBJ)/oi
