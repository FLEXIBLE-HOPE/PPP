! Correction for P1-C1
! Jianghui Geng
! Mar 6 2012
SUBROUTINE corr_p1c1(mjd,sod,cprn,obs,fob)
!!
!*
USE const
IMPLICIT NONE

!*
! The arguments
!!------------------------
INTEGER(IT) ::  mjd
CHARACTER(LEN_PRN) :: cprn(MAXSAT)
REAL(RL) :: sod
REAL(RL) :: obs(MAXSAT,2*MAXFREQ)
CHARACTER(LEN_OBSTYPE) :: fob(MAXSAT,2*MAXFREQ)

  !*
  ! The local variables
  !!------------------------
  INTEGER(IT) :: mjd_sav, isat, i
  REAL(RL) :: cctime,sod_sav,bias(MAXSAT)
  CHARACTER(LEN_PRN) :: prn

  DATA mjd_sav,sod_sav,cctime,bias /0,0.d0,0.d0,MAXSAT*-9.d9/
  SAVE mjd_sav,sod_sav,cctime,bias

  INTEGER(IT) :: pointer_string

  !*
  ! Start the exectuable code
  !!-------------------------------

  !! check and read p1c1bias.2000p
  IF (cctime.EQ.0.d0 .OR.(mjd-cctime)*86400.d0+sod.GT.1800.d0) THEN
    CALL pth_cc2nc(mjd_sav,sod_sav,mjd,sod,bias)
    cctime=mjd+sod/86400.d0
  END IF

  !
  !! correct C1
  CALL pth_cc2nc_corr(cprn,bias,obs,fob)

  RETURN

END SUBROUTINE
