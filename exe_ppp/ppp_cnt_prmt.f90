!!
!! purpose   : count number of the three types of parameters
!!             process, state and determinated parameters
!!
!*
SUBROUTINE ppp_cnt_prmt(CKF,SIT,SAT,NM)
!!
!*
USE info
USE ckdctrl
USE station
USE satellite
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!---------------------
TYPE(CKDCFG) :: CKF
TYPE(SITE) :: SIT(1:*)
TYPE(SATE) :: SAT(MAXSAT)
TYPE(INFM) :: NM(1:*)

  !*
  ! The local variables
  !!-----------------------------
  INTEGER(IT) :: i,isit,isat,isys

  !*
  ! Start the exectuable code
  !!-----------------------------

  DO isit=1, CKF.nsit

    NM(isit).np=0
    NM(isit).nc=0
    NM(isit).ns=0
    NM(isit).npc=0

    !! only position, troposphere and clock parameters are estimated
    IF (SIT(isit).skd(1:1) .EQ. 'S') THEN
      NM(isit).nc=NM(isit).nc+3
    ELSE IF (SIT(isit).skd(1:1) .EQ. 'K') THEN
      NM(isit).np=NM(isit).np+3
    !! kinematic orbits but only white noise position change is estimated
    ELSE IF (SIT(isit).skd(1:2) .EQ. 'DP') THEN
      NM(isit).np=NM(isit).np+3
    !! kinematic orbits
    ELSE IF (SIT(isit).skd(1:2) .EQ. 'DK') THEN
      NM(isit).np=NM(isit).np+3
    !! reduced-dynamic orbit
    !! the dynamic orbit parameters are estimated epoch-wise
    ELSE IF (SIT(isit).skd(1:2) .EQ. 'DE') THEN
      isat=CKF.nprn+SIT(isit).ileo
      IF (SAT(isat).npar .EQ. 0) CYCLE
      NM(isit).np=NM(isit).np+SAT(isat).npar
    ELSE IF (SIT(isit).skd(1:1).NE.'F' .AND. SIT(isit).skd(1:2).NE.'DF') THEN
      WRITE(ERROR_UNIT,'(A)') '***ERROR(ppp_cnt_prmt): unknown station type '//SIT(isit).skd
      CALL exit(1)
    END IF

    !! receiver clock offsets
    IF (CKF%lisb) THEN
      NM(isit).np=NM(isit).np+1
      NM(isit).nc=NM(isit).nc+CKF.nsys-1
    ELSE
      !! ISB is estimated as process paremeter
      NM(isit).np=NM(isit).np+CKF%nsys
    END IF
    !! IFB parameters for GLONASS
    IF (INDEX(CKF.system,'R').NE.0 .AND.  CKF.iref.EQ.INDEX(SYS,'G')) THEN
      NM(isit).nc=NM(isit).nc-1
      DO isat=1, CKF.nprn
        IF (CKF.cprn(isat)(1:1) .EQ. 'R') THEN
          NM(isit).nc=NM(isit).nc+1
        END IF
      END DO
    ELSE IF (CKF.iref .EQ. INDEX(SYS,'R')) THEN
      NM(isit).nc=NM(isit).nc-1
      DO isat=1, CKF.nprn
        IF (CKF.cprn(isat)(1:1) .EQ. 'R') THEN
          NM(isit).nc=NM(isit).nc+1
        END IF
      END DO
    END IF

    !! atmospheric parameters is process parameters
    ! IF (CKF%ztdmod(1:4).NE.'NONE' .AND. CKF%ztdmod(1:3).NE.'FIX' .AND. SIT(isit)%skd(1:1).NE.'D') THEN
    IF (CKF%ztdmod(1:4).NE.'NONE'.AND. SIT(isit)%skd(1:1).NE.'D') THEN
      NM(isit).np=NM(isit).np+1
    END IF

    !! atmospheric grid for north and east
    IF (CKF.grdmod(1:4).NE.'NONE' .AND. SIT(isit).skd(1:1).NE.'D') THEN
      NM(isit).np=NM(isit).np+2
    END IF

    !! slant ionosphere-delay for each satellite
    IF (CKF.cobs(1:3) .EQ. 'RAW') THEN
      ! NM(isit).np=NM(isit).np+CKF.nprn
      DO isat=1,CKF%nprn
        isys=INDEX(SYS,CKF%cprn(isat)(1:1))
        IF (CKF%nfreq(isys).GE.2) THEN
          NM(isit)%np=NM(isit)%np+1
        END IF
      END DO
    END IF
    
    !! the broadcast ephemeris errors absorbing parameter
    IF (CKF%lbds .EQ. .TRUE.) NM(isit)%np=NM(isit)%np+CKF%nprn

    ! @CMT BY XSY: [RECDCB] The third frequency code observation is a useless contribution because its' low wight
    !! the third frequency code bias, only for multi-frequency ppp is need because the code osb is corrected, in some case, it may be estimated for the same system
    IF (CKF%lrecdcb) THEN
      IF (TRIM(CKF.cobs) .EQ. 'RAW') THEN
        DO isat=1, CKF.nprn
          isys=INDEX(SYS,CKF.cprn(isat)(1:1))
          IF (CKF.nfq(isys) .GE. 3) THEN
            NM(isit).nc=NM(isit).nc+CKF.nfq(isys)-2
            ! NM(isit)%np=NM(isit)%np+CKF%nfq(isys)-2
          END IF
        END DO
      END IF
    END IF

    !! estimate the time synchronization bias for LEO satellites as IFB, similar to GLONASS IFB
    IF (CKF%lleoifb) THEN
      NM(isit)%nc=NM(isit)%nc-1
      DO isat=1, CKF%nprn
        IF (CKF%cprn(isat)(1:1) .EQ. 'L') THEN
          NM(isit)%nc=NM(isit)%nc+1
        END IF
      END DO
    END IF

    NM(isit).npc =NM(isit).np+NM(isit).nc
    NM(isit).imtx=NM(isit).npc+NM(isit).ns

    !! check consistence
    IF (NM(isit).imtx .GT. MAXPARSIT) THEN
      WRITE(ERROR_UNIT,'(A)') '***ERROR(ppp_cnt_prmt): too many parameters'
      CALL exit(1)
    END IF

  END DO

  RETURN

END SUBROUTINE
