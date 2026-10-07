! read clocks & observations
! Jianghui Geng
! March 17 2012

SUBROUTINE fcb_rd_clkobs(CKF,SAT,SIT,OB)
!!
!*
USE ckdctrl
USE station
USE satellite
USE observation
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!----------------------
TYPE(CKDCFG) :: CKF
TYPE(SATE) :: SAT(MAXSAT)
TYPE(SITE) :: SIT(MAXSIT)
TYPE(RNXOBS) :: OB(MAXSIT)

  !*
  ! The local variables
  !!------------------------
  LOGICAL(LG) :: lfcb
  INTEGER(IT) :: i,j,isit,isat
  INTEGER(IT) :: nsit,ierr,ifreq(MAXSAT)
  REAL(RL) :: est(MAXSAT),obs(MAXSAT,2*MAXFREQ,MAXSIT)
  CHARACTER(LEN_SITENAME) :: snam(MAXSIT)
  CHARACTER(LEN_ANTENNA) :: antname, antnumb
  CHARACTER(LEN_PRN) :: cprn(MAXSAT)

  !*
  ! The function called
  !!--------------------------
  INTEGER(IT) :: pointer_string

  !*
  ! Start the exectuable code
  !!---------------------------

  !! initialization
  DO isit=1, CKF.nsit
    OB(isit).obs=0.d0
  END DO
  DO isat=1,CKF.nprn
    SAT(isat).sclock=0.d0
  END DO

  !! read satellite clocks and observations from shared memory
  CALL shm_read_clkobs(CKF.mjd,CKF.sod,nsit,snam,cprn,ifreq,est,obs)

  !! store observations
  DO j=1,nsit

    DO isit=1, CKF.nsit
      IF (SIT(isit).name .EQ. snam(j)) EXIT
    END DO

    !! add a new site and initialize it
    IF (isit .GT. CKF.nsit) THEN
      CKF.nsit=CKF.nsit+1
      IF (CKF.nsit .GT. MAXSIT) THEN
        WRITE(ERROR_UNIT,'(A)') '***ERROR(fcb_rd_clkobs): too many stations'
        CALL exit(1)
      END IF
      SIT(isit).name=snam(j)

      CALL read_siteinfo(SIT(isit),CKF.mjd,0.d0,0.d0,ierr)
      IF (ierr .NE. 0) THEN
        WRITE(ERROR_UNIT,'(A,A4)') '***ERROR(fcb_rd_clkobs): no information for ',SIT(isit).name
        CKF.nsit=CKF.nsit-1
        CALL exit(1)
      ELSE
        CALL xyzblh(SIT(isit).x(1:3),1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,SIT(isit).geod)
        CALL rot_enu2xyz(SIT(isit).geod(1),SIT(isit).geod(2),SIT(isit).rot_l2f)
        CALL oceanload_coef(SIT(isit).geod(1),SIT(isit).geod(2),SIT(isit).olc)
        CALL antnam(SIT(isit).name,SIT(isit).anttyp,antname,ierr)
        antnumb = ''
        CALL get_ant_ipt('SITE',CKF.nfreq,CKF.freq,CKF.mjd+CKF.sod/86400.d0,CKF.mjd+CKF.sod/86400.d0,&
                       antname,antnumb,SIT(isit).iptatx,SIT(isit).enu)
      END IF
    END IF


    OB(isit).jd=CKF.mjd
    OB(isit).tsec=CKF.sod
    DO i=1,MAXSAT
      IF (LEN_TRIM(cprn(i)) .NE. 0) THEN
        isat=pointer_string(CKF.nprn,CKF.cprn,cprn(i))
        IF (isat.NE.0 .AND. COUNT(obs(i,1:2*MAXFREQ,j).NE.0.D0).GE.4) THEN
          OB(isit).obs(isat,1:2*MAXFREQ)=obs(i,1:2*MAXFREQ,j)
        END IF
      END IF
    END DO
  END DO

  !! store clocks
  DO i=1,MAXSAT
    IF (LEN_TRIM(cprn(i)) .NE. 0) THEN
      isat=pointer_string(CKF.nprn,CKF.cprn,cprn(i))
      IF (isat .NE. 0) THEN
        SAT(isat).sclock=est(i)
        !! For GLONASS ambiguity resolution in near furture, although we delete the satellites in get_fcbrt_args
        IF (SAT(isat).ifreq .NE. ifreq(i)) THEN
          SAT(isat).ifreq=ifreq(i)
          j=INDEX(SYS,SAT(isat).cprn(1:1))
          CALL get_freq(SAT(isat), CKF.cobs, CKF.nfreq(j), CKF.freq(:,j))
        END IF
        IF (SAT(isat).iptatx .EQ. 0) THEN
          antnumb = CKF.cprn(isat)
          CALL read_svnav(CKF.mjd+CKF.sod/86400.d0,SAT(isat))

          CALL get_ant_ipt(SAT(isat).cprn,CKF.nfreq,CKF.freq,CKF.mjd+CKF.sod/86400.d0,CKF.mjd+CKF.sod/86400.d0,&
                           SAT(isat).type,antnumb,SAT(isat).iptatx,SAT(isat).xyz)
        END IF
      END IF
    END IF
  END DO

  RETURN


END SUBROUTINE
