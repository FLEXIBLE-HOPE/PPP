!
!! purpose  : epoch difference
!! parameter:
!!    input : CKF  -- filter configuration
!!    output: OB   -- observation struct
!! author   : Geng J
!! created  : Dec. 21, 2007
!
SUBROUTINE ppp_epc_diff(isit,CKF,SIT,OB)
!!
!*
USE const
USE station
USE ckdctrl
USE observation
IMPLICIT NONE

!*
! The arguments
!!-----------------------
INTEGER(IT) :: isit
TYPE(SITE) :: SIT
TYPE(CKDCFG) :: CKF
TYPE(RNXOBS) :: OB

  !*
  ! The local variables
  !!----------------------
  LOGICAL(LG) :: lfirst
  INTEGER(IT) :: isat, npar
  REAL(RL) :: tpob(MAXFREQ,MAXSAT,MAXSIT),tpzd(MAXSAT,MAXSIT),tpvar(MAXFREQ,MAXSAT,MAXSIT)
  REAL(RL) :: tmp1(MAXFREQ),tmp2,tmp3(MAXFREQ)

  DATA lfirst/.TRUE./
  SAVE lfirst,tpob,tpzd,tpvar

  !*
  ! Start the exectuable code
  !!----------------------------

  IF (lfirst .EQ. .TRUE.) THEN
    lfirst = .FALSE.
    tpob  = 0.d0
    tpzd  = 0.d0
  END IF

  npar=2+1

  !! form epoch difference and save omc
  DO isat=1, CKF.nprn
    tmp1(1:MAXFREQ)=OB.omc(isat,1:MAXFREQ)
    tmp2=OB.amat(npar,isat)
    tmp3(1:MAXFREQ)=OB.var(isat,1:MAXFREQ)

    IF (OB.flag(isat,1) .NE. 1) THEN
      IF (OB.omc(isat,1).NE.0.d0 .AND. tpob(1,isat,isit).NE.0.d0) THEN
        !! In the same piece, the parameter is same, so the partial party could be subtracted
        IF (SIT.zcor .EQ. 0.d0) THEN
          OB.amat(npar,isat) = OB.amat(npar,isat)-tpzd(isat,isit)
          OB.omc(isat,1:MAXFREQ)= OB.omc(isat,1:MAXFREQ)-tpob(1:MAXFREQ,isat,isit)
        ELSE
          OB.omc(isat,1:MAXFREQ)= OB.omc(isat,1:MAXFREQ)-tpob(1:MAXFREQ,isat,isit)+tpzd(isat,isit)*SIT.zcor
        END IF

        OB.var(isat,1:MAXFREQ)  = OB.var(isat,1:MAXFREQ)+tpvar(1:MAXFREQ,isat,isit)
      ELSE IF (tpob(1,isat,isit) .EQ. 0.d0) THEN
        OB.omc(isat,1:MAXFREQ)  =0.D0
      END IF
    END IF

    tpob(1:MAXFREQ,isat,isit)=tmp1(1:MAXFREQ)
    tpzd(isat,isit)=tmp2
    tpvar(1:MAXFREQ,isat,isit)=tmp3(1:MAXFREQ)
  END DO
  SIT.zcor=0.d0
  SIT.gcor=0.d0

  RETURN

END SUBROUTINE
