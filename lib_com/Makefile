# Makefile for modulesLib Sep 01 13:09:26 2009

all: comLib


PATH_OBJ = ../OBJ_INTEL
FFLAGS   = -w -O3 -fp-model precise -132 -traceback -module $(PATH_OBJ) -mkl=parallel -qopt-zmm-usage=high
FFLAGS7  = -w -O3 -fp-model precise -132 -traceback -module $(PATH_OBJ) -mkl=parallel -qopt-zmm-usage=high
FTN      = ifort
CFLAGS   = -I ./
CC       = icc


LIB = $(PATH_OBJ)/libcom.a

comLib : $(LIB)

$(LIB) : \
	$(LIB)(ambslv.o) \
        $(LIB)(antnam.o) \
	$(LIB)(betau.o) \
	$(LIB)(bds_s2f.o) \
	$(LIB)(brd2xyz.o) \
	$(LIB)(brdtime.o) \
	$(LIB)(bds_code_cor.o) \
	$(LIB)(chi2.o) \
	$(LIB)(chitst.o) \
	$(LIB)(chos.o) \
	$(LIB)(close_all_files.o) \
	$(LIB)(cross.o) \
        $(LIB)(corr_sdb_gxp.o) \
	$(LIB)(code_smooth.o) \
	$(LIB)(day_of_year.o) \
	$(LIB)(dbg_rd_clk.o) \
	$(LIB)(dot.o) \
	$(LIB)(eclips.o) \
	$(LIB)(everett_coeff.o) \
	$(LIB)(everett_interp_orbit.o) \
	$(LIB)(fcb_wgt_mean.o) \
	$(LIB)(ffun.o) \
	$(LIB)(findkey.o) \
	$(LIB)(fund_arg_nutation.o) \
	$(LIB)(get_file_name.o) \
	$(LIB)(get_freq.o) \
	$(LIB)(get_valid_unit.o) \
	$(LIB)(get_wgt_mean.o) \
	$(LIB)(gls2xyz.o) \
	$(LIB)(glsint.o) \
	$(LIB)(gnss_interpolate_orbit.o) \
	$(LIB)(gpsweek.o) \
	$(LIB)(interp_orbit_post.o) \
	$(LIB)(left_justify_string.o) \
	$(LIB)(leo_att_non.o) \
	$(LIB)(leo_everett_interpolate.o) \
	$(LIB)(leo_interpolate_orbit.o) \
	$(LIB)(lower_string.o) \
	$(LIB)(matinv.o) \
	$(LIB)(matmpy.o) \
	$(LIB)(mean_pole.o) \
	$(LIB)(mercator.o) \
	$(LIB)(mjd2date.o) \
	$(LIB)(mjd2doy.o) \
	$(LIB)(mjd2wksow.o) \
	$(LIB)(modified_julday.o) \
	$(LIB)(normal.o) \
	$(LIB)(oceanload_coef.o) \
	$(LIB)(orbit_everett_interpolate.o) \
	$(LIB)(pchi2.o) \
	$(LIB)(phase_windup.o) \
	$(LIB)(pnormal.o) \
	$(LIB)(pointer_int.o) \
	$(LIB)(pointer_string.o) \
        $(LIB)(pointer_string_xx.o) \
	$(LIB)(prob_resol.o) \
	$(LIB)(proc_xl_dirc.o) \
	$(LIB)(pz902wgs84.o) \
	$(LIB)(pv2oe.o) \
	$(LIB)(q2rotmat.o) \
	$(LIB)(rdatx.o) \
	$(LIB)(rdorbh.o) \
	$(LIB)(rdsp3.o) \
	$(LIB)(read_igserp.o) \
	$(LIB)(read_kin_pos.o) \
        $(LIB)(read_kin_pos_lc.o) \
	$(LIB)(read_log_file.o) \
        $(LIB)(read_recfcb.o) \
	$(LIB)(read_rec_clock.o) \
	$(LIB)(read_rnxnav.o) \
	$(LIB)(read_rnxobs.o) \
	$(LIB)(read_rnxobs_head.o) \
	$(LIB)(read_sat_asc.o) \
	$(LIB)(read_sat_clock.o) \
        $(LIB)(read_sdb_gxp.o) \
	$(LIB)(read_siteinfo.o) \
	$(LIB)(read_svnav.o) \
	$(LIB)(read_ztd.o) \
	$(LIB)(read_amb.o) \
	$(LIB)(read_ion.o) \
	$(LIB)(read_ifcb.o) \
	$(LIB)(rot_enu2xyz.o) \
	$(LIB)(rot_scfix2j2000.o) \
	$(LIB)(run_tim.o) \
	$(LIB)(sign_robust.o) \
	$(LIB)(sod2hms.o) \
	$(LIB)(split_string.o) \
        $(LIB)(shcoef_order_wise.o) \
	$(LIB)(table_linear_interpolate.o) \
	$(LIB)(taiutc.o) \
	$(LIB)(timdif.o) \
	$(LIB)(timinc.o) \
	$(LIB)(tim_tag.o) \
	$(LIB)(unit_vector.o) \
	$(LIB)(upper_string.o) \
	$(LIB)(write_ics.o) \
	$(LIB)(write_saved.o) \
	$(LIB)(write_sp3.o) \
	$(LIB)(xyzblh.o) \
	$(LIB)(yeardoy2monthday.o) \
	$(LIB)(yr2year.o) \
	$(LIB)(get_antenna_corr.o) \
	$(LIB)(galileo_att.o) \
	$(LIB)(beidou3_att.o) \
	$(LIB)(mat2quater.o) \


%.o : %.f90
	$(FTN) -c $(FFLAGS) $<
	ar -rv $(LIB) *.o

%.o : %.f
	$(FTN) -c $(FFLAGS7) $<
	ar -rv $(LIB) *.o

clean:
	rm -f $(LIB) *.o
