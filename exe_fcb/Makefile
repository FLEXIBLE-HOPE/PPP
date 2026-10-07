# Makefile for fcb on SunOS 5.5.1 by unimake at Fri Nov 20 13:09:26 MET 1998

all: fcb


PATH_OBJ = ../OBJ_INTEL
FFLAGS   = -w -O3 -xHost -fp-model precise -132 -traceback -module $(PATH_OBJ) -mkl=parallel -qopt-zmm-usage=high
FFLAGS7  = -w -O3 -xHost -fp-model precise -132 -traceback -module $(PATH_OBJ) -mkl=parallel -qopt-zmm-usage=high
FTN      = ifort
CFLAGS   = -I ./
CC       = icc


LIB       = libfcb.a
DE405     = $(PATH_OBJ)/libde405.a
ITRS2GCRS = $(PATH_OBJ)/libitrs2gcrs.a
SOFA      = $(PATH_OBJ)/libsofa.a
ION       = $(PATH_OBJ)/libion.a
TROPO     = $(PATH_OBJ)/libtropo.a
TIDE      = $(PATH_OBJ)/libtide.a
COMM      = $(PATH_OBJ)/libcomm.a
COM       = $(PATH_OBJ)/libcom.a
MODULES   = $(PATH_OBJ)/libmodules.a

fcb : $(LIB) $(DE405) $(ITRS2GCRS) $(SOFA) $(ION) $(TROPO) $(TIDE) $(COMM) $(COM) $(MODULES)
	$(FTN) $(FFLAGS) fcb.f90 $(LIB) $(DE405) $(ITRS2GCRS) $(SOFA) $(ION) $(TROPO) $(TIDE) $(COMM) $(COM) $(MODULES) -pthread -o $(PATH_OBJ)/fcb
	rm -f *.o $(LIB)

$(LIB) : \
	$(LIB)(corr_p1c1.o) \
        $(LIB)(corr_codeosb.o) \
	$(LIB)(dbg_read_clkobs.o) \
	$(LIB)(fcb_add_newamb.o) \
	$(LIB)(fcb_add_obs.o) \
	$(LIB)(fcb_cnt_prmt.o) \
	$(LIB)(fcb_comp_gbeewl.o) \
	$(LIB)(fcb_comp_gbewl.o) \
	$(LIB)(fcb_comp_gbhewl.o) \
	$(LIB)(fcb_comp_gbwl.o) \
	$(LIB)(fcb_comp_gbln.o) \
	$(LIB)(fcb_comp_gbl1.o) \
	$(LIB)(fcb_comp_gfewl.o) \
	$(LIB)(fcb_comp_gfwl.o) \
	$(LIB)(fcb_comp_gfln.o) \
	$(LIB)(fcb_comp_naro.o) \
	$(LIB)(fcb_comp_ewide.o) \
	$(LIB)(fcb_comp_wide.o) \
	$(LIB)(fcb_dat_screen.o) \
   	$(LIB)(fcb_dia.o) \
	$(LIB)(fcb_del_ambi.o) \
	$(LIB)(fcb_ele_hht.o) \
	$(LIB)(fcb.o) \
	$(LIB)(fcb_gps_model.o) \
	$(LIB)(fcb_init.o) \
	$(LIB)(fcb_msu_update.o) \
	$(LIB)(fcb_partial.o) \
	$(LIB)(fcb_rd_clkobs.o) \
	$(LIB)(fcb_rd_rnxoi.o) \
	$(LIB)(fcb_rd_wide.o) \
	$(LIB)(fcb_rt_prep.o) \
	$(LIB)(fcb_rt_sfqc.o) \
	$(LIB)(fcb_scr_adapt.o) \
	$(LIB)(fcb_scr_adapt_xsy.o) \
	$(LIB)(fcb_scr_identify.o) \
	$(LIB)(fcb_scr_identify_xsy.o) \
	$(LIB)(fcb_scr_senvec.o) \
	$(LIB)(fcb_scr_senvec_xsy.o) \
	$(LIB)(fcb_slv_prmt.o) \
	$(LIB)(fcb_solution.o) \
	$(LIB)(fcb_tim_update.o) \
    $(LIB)(fcb_wt_osb.o) \
	$(LIB)(fcb_wt_naro.o) \
	$(LIB)(fcb_wt_wide.o) \
	$(LIB)(fcb_wt_hewide.o) \
	$(LIB)(fcb_wt_eewide.o) \
	$(LIB)(fcb_wt_ewide.o) \
	$(LIB)(fcb_wt_product.o) \
	$(LIB)(get_fcb_args.o) \


%.o : %.f90
	$(FTN) -c $(FFLAGS) $<
	ar -rv $(LIB) *.o

%.o : %.f
	$(FTN) -c $(FFLAGS7) $<
	ar -rv $(LIB) *.o

clean:
	rm -f $(LIB) $(PATH_OBJ)/fcb
