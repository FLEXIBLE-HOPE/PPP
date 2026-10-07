!*
MODULE crd
!!
!!  Consolidated laser Ranging Data format (CRD)
!!      record and variable definitions for FORTRAN
!!      R. Ricklefs UT/CSR July 2007
!!  History:
!!  08/xx/07 - added H3 Target type (v0.27)
!!  11/26/07 - added H4 data quality alert
!!             and #10 stop number
!!             and #20 origin of values (v0.27) rlr.
!!  12/12/07 - commons were re-ordered double precision, integer, character
!!            for efficiency and to keep the compiler happy. :-(
!!  05/07/08 - Expand configuration and data record character fields to
!!             allow up to 40 characters.
!!             - Added detector channel to normalpoint (11) and calibration (40)
!!             records.
!!             - Added field for 'crd' literal to 'h1'.
!!             - Record '21' sky_clarity is now double rather than int.
!!               (v1.00 rlr)
!!
!!
!! ----------------------------------------------------------------------
!!
!!  Header Records
!!
!! H1 - format header
!*
USE par
IMPLICIT NONE


  COMMON /h1/ crd_literal, format_version, &
                 prod_year, prod_mon, prod_day, prod_hour
  CHARACTER*3 :: crd_literal
  INTEGER(IT) :: format_version, prod_year, &
                 prod_mon, prod_day, prod_hour

  ! H2 - station header
  COMMON /h2/ cdp_pad_id, cdp_sys_num, cdp_occ_num, &
                 stn_timescale, stn_name
  CHARACTER*10 :: stn_name
  INTEGER(IT) :: cdp_pad_id, cdp_sys_num, cdp_occ_num, stn_timescale

  ! H3 - spacecraft header
  COMMON /h3/ ilrs_id, sic, norad, SC_timescale, target_type, &
                 target_name
  CHARACTER*10 :: target_name
  INTEGER(IT) :: ilrs_id, sic, norad, SC_timescale, target_type

  ! H4 - Session header
  COMMON /h4/ data_type, start_year, start_mon, start_day, &
              start_hour, start_min, start_sec, &
              end_year, end_mon, end_day, &
              end_hour, end_min, end_sec, &
              data_release, refraction_app_ind, CofM_app_ind, &
              xcv_amp_app_ind, stn_sysdelay_app_ind, &
              SC_sysdelay_app_ind, range_type_ind, &
              data_qual_alert_ind
  INTEGER(IT) :: data_type, start_year, start_mon, start_day
  INTEGER(IT) :: start_hour, start_min, start_sec
  INTEGER(IT) :: end_year, end_mon, end_day, end_hour, end_min, end_sec
  INTEGER(IT) :: data_release, refraction_app_ind, CofM_app_ind
  INTEGER(IT) :: xcv_amp_app_ind, stn_sysdelay_app_ind
  INTEGER(IT) :: SC_sysdelay_app_ind, range_type_ind, data_qual_alert_ind

  ! H8 - End of Session footer

  ! H9 - End of File footer

  !
  !  Configuration Records
  !
  ! C0 - System Configuration Record
  COMMON /c0/ xmit_wavelength, c0_detail_type, config_ids
  INTEGER(IT) :: c0_detail_type
  REAL(RL)    :: xmit_wavelength
  CHARACTER*40 :: config_ids(5)    ! May be more later

  ! C1 - Laser Configuration Record
  COMMON /c1/ prim_wavelength, nom_fire_rate, pulse_energy, &
                 pulse_width, beam_div, &
                 c1_detail_type, pulses_in_semitrain, &
                 laser_type, laser_config_id
  INTEGER(IT) :: c1_detail_type, pulses_in_semitrain
  CHARACTER   :: laser_config_id*40, laser_type*40
  REAL(RL)    :: prim_wavelength, nom_fire_rate, pulse_energy
  REAL(RL)    :: pulse_width, beam_div

  ! C2 - Detector Configuration Record
  COMMON /c2/ app_wavelength, qe, voltage, dark_count, &
                 output_pulse_width, spectral_filter, &
                 spectral_filter_xmission, spatial_filter, &
                 c2_detail_type, signal_proc, &
                 detector_config_id, detector_type, output_pulse_type
  INTEGER(IT) :: c2_detail_type
  CHARACTER   :: detector_config_id*40, detector_type*40
  CHARACTER   :: output_pulse_type*40, signal_proc*40
  REAL(RL)    :: app_wavelength, qe, voltage, dark_count
  REAL(RL)    :: output_pulse_width, spectral_filter
  REAL(RL)    :: spectral_filter_xmission, spatial_filter

  ! C3 - Timing Configuration Record
  COMMON /c3/ c3_detail_type, timing_config_id, time_source, &
                 freq_source, timer, timer_serial_num, epoch_delay_corr
  INTEGER(IT) :: c3_detail_type
  CHARACTER   :: timing_config_id*40, time_source*40
  CHARACTER   :: freq_source*40, timer*40, timer_serial_num*40
  REAL(RL)    :: epoch_delay_corr

  ! C4 - Transponder Configuration Record
  COMMON /c4/ est_stn_utc_offset, est_stn_osc_drift, &
                 est_xponder_utc_offset, est_xponder_osc_drift, &
                 xponder_clock_ref_time, &
                 c4_detail_type, stn_off_drift_app_ind, &
                 SC_off_drift_app_ind, SC_time_simplified_ind, &
                 xponder_config_id
  INTEGER(IT) :: c4_detail_type
  INTEGER(IT) :: stn_off_drift_app_ind, SC_off_drift_app_ind
  INTEGER(IT) :: SC_time_simplified_ind
  CHARACTER   :: xponder_config_id*40
  REAL(RL)    :: est_stn_osc_drift, est_xponder_osc_drift
  ! The next 3 variables sould be "long double" or "real*16"
  REAL(RL)    :: est_stn_utc_offset, est_xponder_utc_offset
  REAL(RL)    :: xponder_clock_ref_time

  !
  ! Data Records
  !
  ! 10 - Range Record
  COMMON /d10/ d10_sec_of_day, d10_time_of_flight, &
                  d10_epoch_event, filter_flag, &
                  d10_detector_channel, stop_number, xcv_amp, &
                  d10_sysconfig_id
  CHARACTER    :: d10_sysconfig_id*40
  REAL(RL)     :: d10_sec_of_day, d10_time_of_flight
  INTEGER(IT)  :: d10_epoch_event, filter_flag, d10_detector_channel
  INTEGER(IT)  :: xcv_amp
  INTEGER(IT)  :: stop_number

  ! 11 - Normal Point Record
  COMMON /d11/ d11_sec_of_day, d11_time_of_flight, &
                  d11_epoch_event, np_window_length, &
                  bin_rms, bin_skew, bin_kurtosis, bin_PmM, &
                  return_rate, d11_sysconfig_id, num_ranges, &
                  d11_detector_channel
  CHARACTER*40 :: d11_sysconfig_id
  REAL(RL)     :: d11_sec_of_day, d11_time_of_flight
  REAL(RL)     :: np_window_length
  REAL(RL)     :: bin_rms, bin_skew, bin_kurtosis, bin_PmM
  REAL(RL)     :: return_rate
  INTEGER(IT)  :: d11_epoch_event, num_ranges
  INTEGER(IT)  :: d11_detector_channel

  ! 12 - Range Supplement Record
  COMMON /d12/ d12_sec_of_day, refraction_corr, &
                  target_CofM_corr, nd_value, time_bias, &
                  d12_sysconfig_id
  CHARACTER*40 :: d12_sysconfig_id
  REAL(RL)     :: d12_sec_of_day, refraction_corr
  REAL(RL)     :: target_CofM_corr, nd_value, time_bias

  ! 20 - Meteorological Record
  COMMON /d20/ d20_sec_of_day, pressure, &
                  temperature, humidity, value_origin
  REAL(RL)     :: d20_sec_of_day, pressure, temperature, humidity
  INTEGER(IT)  :: value_origin

  ! 21 - Meteorological Supplement Record
  COMMON /D21/ d21_sec_of_day, wind_speed, wind_direction, &
                  sky_clarity, visibility, atmospheric_seeing, &
                  cloud_cover, precip_type
  CHARACTER    :: precip_type*40
  REAL(RL)     :: d21_sec_of_day, wind_speed, wind_direction
  REAL(RL)     :: sky_clarity
  INTEGER(IT)  :: visibility, atmospheric_seeing, cloud_cover

  ! 30 - Pointing Angles Record
  COMMON /D30/ d30_sec_of_day, azimuth, elevation, direction_ind, &
                  angle_origin_ind, refraction_corr_ind
  REAL(RL)     :: d30_sec_of_day, azimuth, elevation
  INTEGER(IT)  :: direction_ind, angle_origin_ind, refraction_corr_ind

  ! 40 - Calibration Record
  COMMON /d40/ d40_sec_of_day, type_of_data, d40_sysconfig_id, &
                  num_points_recorded, num_point_used, &
                  one_way_target_dist, cal_sys_delay, cal_delay_shift, &
                  cal_rms, cal_skew, cal_kurtosis, cal_PmM, &
                  cal_type_ind, cal_shift_type_ind, &
                  d40_detector_channel
  CHARACTER*40 :: d40_sysconfig_id
  REAL(RL)     :: d40_sec_of_day
  REAL(RL)     :: one_way_target_dist, cal_sys_delay, cal_delay_shift
  REAL(RL)     :: cal_rms, cal_skew, cal_kurtosis, cal_PmM
  INTEGER(IT)  :: type_of_data, num_points_recorded, num_point_used
  INTEGER(IT)  :: cal_type_ind, cal_shift_type_ind, d40_detector_channel

  ! 50 - Session Statistics Record
  COMMON /d50/ sess_rms, sess_skew, sess_kurtosis, sess_PmM, &
                  data_qual_ind, d50_sysconfig_id
  CHARACTER*40 :: d50_sysconfig_id
  REAL(RL)     :: sess_rms, sess_skew, sess_kurtosis, sess_PmM
  INTEGER(IT)  :: data_qual_ind

  ! 60 - Compatibility Record
  COMMON /d60/ sys_change_ind,sys_config_ind, d60_sysconfig_id
  CHARACTER*40 :: d60_sysconfig_id
  INTEGER(IT)  :: sys_change_ind,sys_config_ind

  ! 9X - User Defined Record

  ! 00 - Comment Record
  COMMON /d00/ comment
  CHARACTER*80 comment

END MODULE
