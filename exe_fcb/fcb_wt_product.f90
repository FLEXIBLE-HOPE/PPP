! Write satellite products for PPP
! Jianghui Geng
! March 19 2012

SUBROUTINE fcb_wt_product(CKF,SAT,UPD)
!!
!*
USE ckdctrl
USE satellite
USE ambiguity
IMPLICIT NONE

!*
! The arguments
!!----------------------
TYPE(SATE) :: SAT(MAXSAT)
TYPE(CKDCFG) :: CKF
TYPE(FCB) :: UPD

  !*
  ! The local variables
  !!------------------------
  INTEGER(IT) :: isat,isit,ifreq(MAXSAT)
  CHARACTER(LEN_PRN) :: cprn(MAXSAT)
  REAL(RL) :: est(MAXSAT)

  !*
  ! Start the exectuable code
  !!--------------------------

  cprn=''
  est=0.d0
  DO isat=1,CKF.nprn
    ifreq(isat)=SAT(isat).ifreq
    cprn(isat)=CKF.cprn(isat)
    est(isat)=SAT(isat).sclock
  END DO
  CALL shm_writer(CKF.mjd,CKF.sod,cprn,ifreq,est,UPD.wfcb,UPD.wsl,UPD.nfcb,UPD.nsl)

  RETURN

END SUBROUTINE
