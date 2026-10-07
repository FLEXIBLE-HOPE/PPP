# Makefile for ppp on SunOS 5.5.1 by unimake at Fri Nov 20 13:09:26 MET 1998

all: ppp


PATH_OBJ = ../OBJ_INTEL
FFLAGS   = -w -O3 -132 -traceback -mkl  -module $(PATH_OBJ)
FFLAGS7  = -w -O3 -132 -traceback -mkl  -module $(PATH_OBJ)
FTN      = ifort
CFLAGS   = -I ./
CC       = icc


LIB       = libppp.a
OI        = $(PATH_OBJ)/liboi.a
DE405     = $(PATH_OBJ)/libde405.a
ITRS2GCRS = $(PATH_OBJ)/libitrs2gcrs.a
SOFA      = $(PATH_OBJ)/libsofa.a
ION       = $(PATH_OBJ)/libion.a
TROPO     = $(PATH_OBJ)/libtropo.a
TIDE      = $(PATH_OBJ)/libtide.a
COMM      = $(PATH_OBJ)/libcomm.a
COM       = $(PATH_OBJ)/libcom.a
LAMBDA    = $(PATH_OBJ)/liblambda.a
MODULES   = $(PATH_OBJ)/libmodules.a

MKLROOT = /opt/intel/oneapi/mkl/2024.2/
Link_mul = ${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
ComFLAG = -i8 -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include


ppp : $(LIB) $(OI) $(DE405) $(ITRS2GCRS) $(SOFA) $(ION) $(TROPO) $(TIDE) $(COMM) $(COM) $(LAMBDA) $(MODULES)
	$(FTN) $(FFLAGS) ppp.f90 $(LIB) $(OI) $(DE405) $(ITRS2GCRS) $(SOFA) $(ION) $(TROPO) $(TIDE) $(COMM) $(COM) $(LAMBDA) $(MODULES) ${Link_mul} -pthread -o $(PATH_OBJ)/ppp
	rm -f *.o $(LIB)

$(LIB) : \
	$(LIB)(chos_sort.o) \
	$(LIB)(chsip.o) \
	$(LIB)(corr_hisi.o) \
	$(LIB)(corr_iondel_grid.o) \
	$(LIB)(corr_iondel_sd.o) \
	$(LIB)(corr_osb_code.o) \
	$(LIB)(corr_osb.o) \
	$(LIB)(corr_osb_mosb.o) \
	$(LIB)(corr_osb_upd.o) \
	$(LIB)(corr_p1c1.o) \
	$(LIB)(cyslip.o) \
	$(LIB)(dbg_rd_fcb.o) \
	$(LIB)(det.o) \
	$(LIB)(dop.o) \
	$(LIB)(get_ppp_args.o) \
	$(LIB)(ppp_abfix_emw.o) \
	$(LIB)(ppp_abfix_fcb_pri.o) \
	$(LIB)(ppp_abfix_fcb_uduc_iar.o) \
	$(LIB)(ppp_abfix_fcb_uduc_war.o) \
	$(LIB)(ppp_abfix_gfewl.o) \
	$(LIB)(ppp_abfix_gfwl.o) \
	$(LIB)(ppp_abfix_l1.o) \
	$(LIB)(ppp_abfix_lambda.o) \
	$(LIB)(ppp_abfix_mf.o) \
	$(LIB)(ppp_abfix_mw.o) \
	$(LIB)(ppp_abfix_pri_dualsd.o) \
	$(LIB)(ppp_abfix_pri.o) \
	$(LIB)(ppp_abfix_uduc_iar.o) \
	$(LIB)(ppp_abfix_uduc_war.o) \
	$(LIB)(ppp_abfix_wlpar.o) \
	$(LIB)(ppp_add_ambcon.o) \
	$(LIB)(ppp_add_newamb.o) \
	$(LIB)(ppp_addob_pri.o) \
	$(LIB)(ppp_amb_fcb_uduc_sd.o) \
	$(LIB)(ppp_amb_uduc_sd.o) \
	$(LIB)(ppp_cmpt_lglw.o) \
	$(LIB)(ppp_cnt_prmt.o) \
	$(LIB)(ppp_con_bias.o) \
	$(LIB)(ppp_dat_screen.o) \
	$(LIB)(ppp_del_ambi.o) \
	$(LIB)(ppp_del_par.o) \
	$(LIB)(ppp_det_cycle.o) \
	$(LIB)(ppp_dia.o) \
	$(LIB)(ppp_ele_hht.o) \
	$(LIB)(ppp_epc_diff.o) \
	$(LIB)(ppp.o) \
	$(LIB)(ppp_fright.o) \
	$(LIB)(ppp_gps_model.o) \
	$(LIB)(ppp_igg3.o) \
	$(LIB)(ppp_init_pri.o) \
	$(LIB)(ppp_map_invx.o) \
	$(LIB)(ppp_map_raw.o) \
	$(LIB)(ppp_msu_update.o) \
	$(LIB)(ppp_oi.o) \
	$(LIB)(ppp_out_omc.o) \
	$(LIB)(ppp_partial.o) \
	$(LIB)(ppp_plc_fixed.o) \
	$(LIB)(ppp_plc_float.o) \
	$(LIB)(ppp_plc_l1.o) \
	$(LIB)(ppp_rd_clkfcb.o) \
	$(LIB)(ppp_rd_rnxoi.o) \
	$(LIB)(ppp_rkf.o) \
	$(LIB)(ppp_rt_prep.o) \
	$(LIB)(ppp_scr_adapt.o) \
	$(LIB)(ppp_scr_grid.o) \
	$(LIB)(ppp_scr_identify.o) \
	$(LIB)(ppp_scr_senvec.o) \
	$(LIB)(ppp_slv_inform.o) \
	$(LIB)(ppp_solution.o) \
	$(LIB)(ppp_spp.o) \
	$(LIB)(ppp_tim_update_block.o) \
	$(LIB)(ppp_tim_update.o) \
	$(LIB)(ppp_wt_cfg.o) \
	$(LIB)(read_osb.o) \
	$(LIB)(read_osb_mosb.o) \
	$(LIB)(read_osb_upd.o) \
	$(LIB)(rep_integer.o) \
	$(LIB)(rep_l1phase.o) \
	$(LIB)(rep_rawphase.o) \
	$(LIB)(rep_widelane.o) \

%.o : %.f90
	$(FTN) -c ${ComFLAG} $(FFLAGS) $<
	ar -rv $(LIB) *.o
 
%.o : %.f
	$(FTN) -c ${ComFLAG} $(FFLAGS7) $<
	ar -rv $(LIB) *.o

clean:
	rm -f $(LIB) $(PATH_OBJ)/ppp
