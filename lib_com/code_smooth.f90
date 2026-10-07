!*
SUBROUTINE phase_smooth_code(CKF,OB,SAT,isit,win)
!!
!! Hatch filter
!!
!*
USE const
USE ckdctrl
USE satellite
USE observation
IMPLICIT NONE

!*
! The arguments
!!---------------------------
INTEGER(IT) :: win,isit
TYPE(RNXOBS) :: OB
TYPE(SATE) :: SAT(MAXSAT)
TYPE(CKDCFG) :: CKF


  !*
  ! The local variables
  !!---------------------------
  INTEGER(IT) :: isat,ifreq,isys,num,i
  INTEGER(IT) :: c(MAXSIT,MAXSAT,MAXFREQ)=0
  ! mean bias of phase and code
  REAL(RL) :: sobs(MAXSIT,MAXSAT,MAXFREQ)=0.d0

  LOGICAL(LG) :: lnarc
  REAL(RL) :: pcb,avg

  SAVE c,sobs

  !*
  ! Start the exectuable code
  !!---------------------------

  DO isat=1, CKF.nprn

    isys=INDEX(SYS,CKF.cprn(isat)(1:1))

    DO ifreq=1, CKF.nfreq(isys)

      ! phase-smoothed code, code and phase are exist
      IF (OB.obs(isat,ifreq).NE.0.d0 .AND. OB.obs(isat,ifreq+MAXFREQ).NE.0.d0) THEN

        lnarc=.FALSE.
        IF (CKF.nfreq(isys) .NE. 2) THEN
          pcb=OB.obs(isat,ifreq+MAXFREQ)-OB.obs(isat,ifreq)*VEL_LIGHT/SAT(isat).freq(ifreq)
          IF (OB.flag(isat,ifreq) .EQ. 1) lnarc=.TRUE.
        ELSE
          IF (ifreq .EQ. 1) THEN
            pcb=OB.obs(isat,ifreq+MAXFREQ)-OB.obs(isat,ifreq)*VEL_LIGHT/SAT(isat).freq(ifreq)- &
                2.d0*SAT(isat).fac(2)*(OB.obs(isat,1)*VEL_LIGHT/SAT(isat).freq(1)-OB.obs(isat,2)*VEL_LIGHT/SAT(isat).freq(2))
          ELSE IF (ifreq .EQ. 2) THEN
            pcb=OB.obs(isat,ifreq+MAXFREQ)-OB.obs(isat,ifreq)*VEL_LIGHT/SAT(isat).freq(ifreq)- &
                2.d0*SAT(isat).fac(1)*(OB.obs(isat,1)*VEL_LIGHT/SAT(isat).freq(1)-OB.obs(isat,2)*VEL_LIGHT/SAT(isat).freq(2))
          END IF
          IF (OB.flag(isat,ifreq).EQ.1 .OR. OB.flag(isat,ifreq+1).EQ.1) lnarc=.TRUE.
        END IF


        ! new ambiguity arc, we should restart the smooth
        IF (lnarc .EQ. .TRUE.) THEN
     !   IF (OB.flag(isat,ifreq) .EQ. 1) THEN
          ! valied epochs in the arc/win
          c(isit,isat,ifreq)=1
          sobs(isit,isat,ifreq)=pcb
        ! old ambiguity arc
        ELSE
          c(isit,isat,ifreq)=c(isit,isat,ifreq)+1
          num=c(isit,isat,ifreq)
          avg=sobs(isit,isat,ifreq)
          ! mean value
          IF (num .GE. win) THEN
            avg=1.d0/win*pcb+(win-1.d0)/win*avg
          ELSE
            avg=1.d0/num*pcb+(num-1.d0)/num*avg
          END IF
          ! save the upgrated average
          sobs(isit,isat,ifreq)=avg
        END IF

        OB.obs(isat,ifreq+MAXFREQ)=OB.obs(isat,ifreq)*VEL_LIGHT/SAT(isat).freq(ifreq)+sobs(isit,isat,ifreq)

      ! as we already flag the positioning of ambiguity, hence, for the bad observations, we do not restart the smooth, just keep it
      END IF
      
    END DO

  END DO

  RETURN

END SUBROUTINE
