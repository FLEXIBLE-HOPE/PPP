!
!! purpose  : read an observation from BNC2.9
!! parameter:
!!    input : SAT -- satellites
!!    output: SIT -- sites
!! author   : Jianghui Geng
!! created  : Feb 14 2012
!! update   : June 2013
!
SUBROUTINE fcb_rd_rnxoi(CKF,SAT,SIT,OB)
!!
!*
USE par
USE station
USE ckdctrl
USE satellite
USE observation
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!------------------------
TYPE(CKDCFG) :: CKF
TYPE(SATE) :: SAT(MAXSAT)
TYPE(SITE) :: SIT(MAXSIT)
TYPE(RNXOBS) :: OB(MAXSIT)

  !*
  ! The local variables
  !!-----------------------
  INTEGER(IT) :: i = 0
  INTEGER(IT) :: k = 0
  INTEGER(IT) :: ic = 0
  INTEGER(IT) :: jc = 0
  INTEGER(IT) :: kc = 0
  INTEGER(IT) :: ilen = 0
  INTEGER(IT) :: ird = 0
  INTEGER(IT) :: ierr = 0
  INTEGER(IT) :: mjd = 0
  INTEGER(IT) :: week = 0
  INTEGER(IT) :: isit = 0
  INTEGER(IT) :: isat = 0
  INTEGER(IT) :: fepoch = 0
  INTEGER(IT) :: ifreq = 0
  INTEGER(IT) :: isys = 0
  INTEGER(IT) :: iobs = 0

  REAL(RL) :: sod = 0.D0
  REAL(RL) :: wks = 0.D0

  CHARACTER(LEN_ANTENNA) :: antname = ''
  CHARACTER(LEN_ANTENNA) :: antnumb = ''
  CHARACTER(LEN=1024) :: rdone = ' '
  CHARACTER(LEN=2048) :: rdres = ' '
  SAVE ilen, rdres

  CHARACTER(LEN_OBSTYPE) :: code

  !*
  ! The function called
  !!------------------------
  INTEGER(IT) :: pointer_string

  !*
  ! Start the exectuable code
  !!--------------------------

  !! zero obs
  DO isit=1, CKF.nsit
    OB(isit).obs = 0.d0
    OB(isit).fob = ' '
    OB(isit).nprn = CKF.nprn
    OB(isit).cprn = CKF.cprn
  END DO

  CKF.mjd = 0
  CKF.sod = 0.d0

  !! read a full epoch of measurements
  fepoch=0
  DO WHILE(fepoch .EQ. 0)

    !! read if there is no one full record
    i = INDEX(rdres(1:ilen),'%')
    IF (i.EQ.0 .OR. i.EQ.ilen) THEN

      !! get data from opened socket
      IF (ilen+1024 .GT. 2048) THEN
        WRITE(ERROR_UNIT,'(A,I5)') '***ERROR(fcb_rd_rnxoi): too short character array ', ilen
        WRITE(ERROR_UNIT,'(A)') rdres(1:ilen)
        CALL exit(1)
      END IF

      IF (CKF.ldebug .EQ. .TRUE.) THEN
        CALL read_saved(rdone,ird)
      ELSE
        CALL read_rtcm3(CKF.sock,rdone,CKF.port,ird)
        !!CALL write_saved(rdone,ird)
      END IF

      rdres=rdres(1:ilen)//rdone(1:ird)
      ilen=ilen+ird
    END IF

    !! grab measurements
    kc=1
    DO WHILE(.TRUE.)
      ic=INDEX(rdres(kc:ilen),'%')
      !! new epoch
      jc=INDEX(rdres(kc:ilen),'%%')

      IF (ic .NE. jc) jc=0
      IF (ic+kc-1 .EQ. ilen) ic=0  ! leave % in case of %%

      IF (ic.NE.0 .AND. fepoch.EQ.0) THEN

        ic=ic+kc-1

        !! new epoch
        IF (rdres(kc:kc) .EQ. '>') THEN
          !READ(rdres(kc+2:),'(I4,F15.7)',err=51) week, wks !old bnc version
          READ(rdres(kc+1:ic-1),*,err=51) week, wks
          IF (wks.LT.0.d0 .OR. wks.GT.604800.d0) goto 51
          mjd=week*7+44243+INT(wks/86400.d0)+1
          sod=dmod(wks,86400.d0)
          IF (CKF.mjd .EQ. 0) THEN
            CKF.mjd=mjd
            CKF.sod=sod
          ELSE IF(CKF.mjd.NE.mjd .OR. CKF.sod.NE.sod) THEN
            fepoch=1
            CYCLE
          END IF
          GOTO 51
        END IF

        !! search existing sites
        DO isit=1,CKF.nsit
          IF (SIT(isit).name .EQ. rdres(kc:kc+3)) EXIT
        END DO

        !! add a new site and initialize it
        IF (isit .GT.CKF.nsit) THEN
          CKF.nsit=CKF.nsit+1
          IF (CKF.nsit .GT. MAXSIT) THEN
            WRITE(ERROR_UNIT,'(A)') '***ERROR(ckd_rd_rnxoi): too many stations'
            CALL exit(1)
          ENDIF
          SIT(isit).name=rdres(kc:kc+3)

          CALL read_siteinfo(SIT(isit),CKF.mjd,0.d0,0.d0,ierr)
          IF (ierr .NE. 0) THEN
            WRITE(OUTPUT_UNIT,'(A,A4)') '%%%MESSAGE(ckd_rd_rnxoi): no information for ',SIT(isit).name
            CKF.nsit=CKF.nsit-1
            GOTO 51
          ELSE
            CALL xyzblh(SIT(isit).x(1:3),1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,SIT(isit).geod)
            CALL rot_enu2xyz(SIT(isit).geod(1),SIT(isit).geod(2),SIT(isit).rot_l2f)
            CALL oceanload_coef(SIT(isit).geod(1),SIT(isit).geod(2),SIT(isit).olc)
            CALL antnam(SIT(isit).name,SIT(isit).anttyp,antname,ierr)
            antnumb = ''
            CALL get_ant_ipt('SITE',CKF.nfreq,CKF.freq,CKF.mjd+CKF.sod/86400.d0,CKF.mjd+CKF.sod/86400.d0,&
                           antname,antnumb,SIT(isit).iptatx,SIT(isit).enu)
          END IF
          OB(isit).obs = 0.D0
          OB(isit).fob = ' '
        END IF
        OB(isit).jd=CKF.mjd
        OB(isit).tsec=CKF.sod

        !! recognize satellite
        i =INDEX(rdres(kc:ic-1),' ')+kc-1
        IF (INDEX(CKF.system,rdres(i+1:i+1)) .EQ. 0) GOTO 51

        isat=pointer_string(CKF.nprn,CKF.cprn,rdres(i+1:i+3))
        !! whether satellite is incorporated
        IF (isat .NE. 0) THEN

          isys = INDEX(SYS, rdres(i+1:i+1))

          !IF (rdres(i+1:i+1) .EQ. 'R') THEN
          !  READ(rdres(i+5:),*,ERR=51) ifreq
          !  IF (SAT(isat).ifreq .NE. ifreq) THEN
          !    SAT(isat).ifreq = ifreq
          !    CALL get_freq(SAT(isat), CKF.cobs, CKF.nfreq(isys), CKF.freq(:,isys))
          !  END IF
          !END IF

          DO ifreq=1, CKF.nfreq(isys)
            DO iobs=1, LEN(OBSTYPE)
              WRITE(code,'(A1,A1,A1)') 'C', CKF.freq(ifreq,isys)(2:2),OBSTYPE(iobs:iobs)
              i =INDEX(rdres(kc:ic-1),code)
              IF (i .EQ. 0) CYCLE
              i = i+kc-1
              READ(rdres(i+3:),*,iostat=ierr) OB(isit).obs(isat,MAXFREQ+ifreq)
              IF (ierr .NE. 0) then
                OB(isit).obs(isat,MAXFREQ+ifreq) = 0.d0
                OB(isit).fob(isat,MAXFREQ+ifreq) = ' '
                CYCLE
              ELSE
                OB(isit).fob(isat,MAXFREQ+ifreq) = code
              END IF

              WRITE(code,'(A1,A1,A1)') 'L', CKF.freq(ifreq,isys)(2:2),OBSTYPE(iobs:iobs)
              i =INDEX(rdres(kc:ic-1),code)
              IF (i .EQ. 0) CYCLE
              i = i+kc-1
              READ(rdres(i+3:),*,iostat=ierr) OB(isit).obs(isat,ifreq)
              IF (ierr .NE. 0) then
                OB(isit).obs(isat,ifreq) = 0.d0
                OB(isit).fob(isat,ifreq) = ' '
                CYCLE
              ELSE
                OB(isit).fob(isat,ifreq) = code
              END IF

              IF (OB(isit).obs(isat,MAXFREQ+ifreq).NE.0.D0 .AND. OB(isit).obs(isat,ifreq).NE.0.D0) EXIT
            END DO
          END DO
          IF (COUNT(OB(isit).obs(isat,1:2*MAXFREQ).NE.0.D0) .LT. 4) THEN
            OB(isit).obs(isat,1:2*MAXFREQ) = 0.D0
            OB(isit).fob(isat,1:2*MAXFREQ) = ' '
          END IF

          IF (SAT(isat).iptatx .EQ. 0) THEN
            antnumb = CKF.cprn(isat)
            CALL read_svnav(CKF.mjd+CKF.sod/86400.d0,SAT(isat))

            CALL get_freq(SAT(isat), CKF.cobs, CKF.nfreq(isys), CKF.freq(:,isys))

            CALL get_ant_ipt(SAT(isat).cprn,CKF.nfreq,CKF.freq,CKF.mjd+CKF.sod/86400.d0,CKF.mjd+CKF.sod/86400.d0,&
                             SAT(isat).type,antnumb,SAT(isat).iptatx,SAT(isat).xyz)
          END IF

        END IF

        !! move pointer to the next measurement
51      CONTINUE

        !! a full epoch has been read, should return now
        IF (jc .NE. 0) THEN
          kc=ic+2
          fepoch=1
        ELSE
          kc=ic+1
        END IF

      ELSE
        !! move the remaining string to the head
        rdres(1:ilen-kc+1)=rdres(kc:ilen)
        ilen=ilen-kc+1
        rdres(ilen+1:2048)=' '
        EXIT
      ENDIF

    !! next measurement
    END DO

  !! next port read
  END DO

  !! correction for P1-C1
  DO isit=1,CKF.nsit
    CALL corr_p1c1(CKF.mjd,CKF.sod,CKF.cprn,OB(isit).obs,OB(isit).fob)
  END DO

  RETURN

END SUBROUTINE
