!*
!! purpose  : detect bad data/cycle slips for single frequency
!!
!
SUBROUTINE fcb_rt_sfqc(CKF,OB,SAT)
!!
!*
USE const
USE ckdctrl
USE satellite
USE observation
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! Start the exectuable
!!------------------------
TYPE(CKDCFG) :: CKF
TYPE(RNXOBS) :: OB(1:*)
TYPE(SATE) :: SAT(MAXSAT)

  !*
  ! The local variables
  !!------------------------
  LOGICAL(LG) :: lfirst
  INTEGER(IT) :: isit,isat,isys,ifreq,ljd(MAXFREQ,MAXSAT,MAXSIT)
  REAL(RL) :: ltm(MAXFREQ,MAXSAT,MAXSIT),lsg(MAXFREQ,MAXSAT,MAXSIT)
  REAL(RL) :: dt,lgr

  DATA lfirst /.TRUE./
  SAVE lfirst,ljd,ltm,lsg

  !*
  ! The function called
  !!------------------------
  REAL(RL) :: timdif


  !*
  ! Start the exectuable code
  !!----------------------------

  IF (lfirst) THEN
    lfirst=.FALSE.
    ljd=0
    ltm=0.d0
    lsg=0.d0
  END IF

  !! real-time preprocessing
  DO isit=1, CKF.nsit
    DO isat=1, CKF.nprn

      isys=INDEX(SYS,CKF.cprn(isat)(1:1))


      DO ifreq=1, CKF.nfq(isys)

        IF (OB(isit).obs(isat,ifreq).EQ.0.d0 .OR. OB(isit).obs(isat,MAXFREQ+ifreq).EQ.0.d0) CYCLE

        !! gap is so large that a new ambiguity has to be set
        IF (ljd(ifreq,isat,isit) .NE. 0) THEN
          dt=timdif(CKF.mjd,CKF.sod,ljd(ifreq,isat,isit),ltm(ifreq,isat,isit))
          IF (dt .GT. CKF.gap) OB(isit).flag(isat,ifreq)=1
        END IF

        !! LG ionosphere observation, widelane ambiguity
        lgr=OB(isit).obs(isat,ifreq)*SAT(isat).lamda(ifreq)-OB(isit).obs(isat,MAXFREQ+ifreq)

        !! new ambiguity
        IF (ljd(ifreq,isat,isit) .NE. 0) THEN
          !IF (DABS(lgr-lsg(ifreq,isat,isit)) .GT. CKF.lg*DSQRT(dt)) THEN
          !  OB(isit).flag(isat,ifreq)=1
          !END IF
        ELSE
          !! first ambiguity for this satellite
          OB(isit).flag(isat,ifreq)=1
        END IF

        !! reserve obs
        ljd(ifreq,isat,isit)=CKF.mjd
        ltm(ifreq,isat,isit)=CKF.sod
        lsg(ifreq,isat,isit)=lgr
      END DO
    END DO
  END DO

  RETURN

END SUBROUTINE
