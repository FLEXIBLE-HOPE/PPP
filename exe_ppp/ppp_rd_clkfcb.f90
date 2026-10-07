! Read satellite clocks and FCBs
! Jianghui Geng
! March 19 2012

SUBROUTINE ppp_rd_clkfcb(mjd,sod,CKF,SAT,UPD)
!*
USE ckdctrl
USE satellite
USE ambiguity
IMPLICIT NONE

!*
! The arguments
!!-------------------------
INTEGER(IT) :: mjd
REAL(RL) :: sod
TYPE(CKDCFG) :: CKF
TYPE(SATE) :: SAT(MAXSAT)
TYPE(FCB) :: UPD

  !*
  ! The local variables
  !!--------------------------
  INTEGER(IT) :: i,isat,ifreq(MAXSAT)
  CHARACTER(LEN_PRN) :: cprn(MAXSAT)
  REAL(RL) :: est(MAXSAT)
  REAL(RL) :: fwl(MAXSAT),swl(MAXSAT),fnl(MAXSAT),snl(MAXSAT)

  INTEGER(IT) :: date_time(8)
  CHARACTER(LEN=12) :: real_clock(3)
  INTEGER(IT) :: mjd0, mjd1
  REAL(RL) :: sod0, sod1

  !*
  ! The function used
  !!---------------------------
  LOGICAL(LG) :: pth_acqcf_read
  INTEGER(IT) :: pointer_string
  INTEGER(IT) :: modified_julday
  REAL(RL) :: timdif

  !*
  ! Start the exectuable code
  !!---------------------------

  !! initialization
  DO isat=1, CKF.nprn
    SAT(isat).sclock=0.d0
  END DO
  UPD.nfcb=10.d0
  UPD.nsl=10.d0

  !! read satellite clocks from shared memory

  CALL DATE_AND_TIME(real_clock(1),real_clock(2),real_clock(3),date_time)
  mjd0=modified_julday(date_time(3),date_time(2),date_time(1))
  sod0=date_time(5)*3600.d0+date_time(6)*60.d0+date_time(7)
  mjd1=mjd0
  sod1=sod0

  !! in the interval, the clock is good for use
  DO WHILE(timdif(mjd1,sod1,mjd0,sod0).LE. CKF.dintv)
    CALL DATE_AND_TIME(real_clock(1),real_clock(2),real_clock(3),date_time)
    mjd1=modified_julday(date_time(3),date_time(2),date_time(1))
    sod1=date_time(5)*3600.d0+date_time(6)*60.d0+date_time(7)

  IF (pth_acqcf_read(mjd,sod,cprn,ifreq,est,fwl,swl,fnl,snl)) THEN

    !! save clocks and fcbs
    DO i=1, MAXSAT
      IF (LEN_TRIM(cprn(i)) .NE. 0) THEN
        isat=pointer_string(CKF.nprn,CKF.cprn,cprn(i))
        IF (isat .NE. 0) THEN
          SAT(isat).sclock=est(i)
          IF (fwl(i) .NE. 10.d0) THEN
            UPD.wfcb(isat)=fwl(i)
            UPD.wsl(isat)=swl(i)
          END IF
          IF (fnl(i).NE.10.d0 .AND. snl(i).LT.0.1d0) THEN
            UPD.nfcb(isat)=fnl(i)
            UPD.nsl(isat)=snl(i)
          END IF
        END IF
      END IF
    END DO

    RETURN

  END IF
  END DO

  RETURN

END SUBROUTINE
