!*
SUBROUTINE read_log_file(iepo,jd0,sod0,dintv0,nprn,cprn,SIT,OB)
!!
!*
USE par
USE station
USE observation
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!------------------------
INTEGER(IT) :: iepo,jd0,nprn
CHARACTER(LEN_PRN) :: cprn(MAXSAT)
REAL(RL) :: sod0,dintv0
TYPE(SITE) :: SIT
TYPE(RNXOBS) :: OB

  !*
  ! Local variable
  !!----------------------------
  INTEGER(IT) :: i,lfn,jepo,kepo
  INTEGER(IT) :: ifreq,id,ierr
  REAL(RL) :: dt,dintv,deljd,ambjd
  LOGICAL(LG) :: lexist
  CHARACTER(3) :: cflag,aprn
  CHARACTER(LEN_STRING) :: line

  !*
  ! The function called
  !!-----------------------------
  INTEGER(IT) :: get_valid_unit
  INTEGER(IT) :: pointer_string
  INTEGER(IT) :: pointer_int

  !*
  ! Start of executable code
  !!------------------------

  !! check and open log file
  IF(SIT.logfile(1:1) .EQ. ' ') RETURN
  IF(SIT.logfile(1:1).NE.' ' .AND. SIT.lfnlog.EQ.0) THEN
    INQUIRE(file=SIT.logfile,exist=lexist)
    IF(.NOT.lexist) THEN
      RETURN
    END IF

    SIT.lfnlog=get_valid_unit(10)
    OPEN(SIT.lfnlog,FILE=SIT.logfile)
    READ(SIT.lfnlog,'(26x,i6,f10.3,f10.3)',iostat=ierr) SIT.jdlog,SIT.sodlog,SIT.intvlog
    IF(ierr.ne.0) GOTO 200
    READ(SIT.lfnlog,'(a)') line
    DO WHILE(line(1:14).ne.'%End of header')
      READ(SIT.lfnlog,'(a)') line
    END DO
    IF(iepo .EQ. 0) RETURN
  END IF

  DO WHILE(.TRUE.)
    READ(SIT.lfnlog,'(a3,x,a3,2i7,i4)',END=100,err=200) cflag,aprn,jepo,kepo,id
    dt=(jd0-SIT.jdlog)*86400.d0+(sod0-SIT.sodlog)+(iepo-1)*dintv0-(jepo-1)*SIT.intvlog
    IF(dt .LT. -0.5*SIT.intvlog) THEN
      BACKSPACE SIT.lfnlog
      GOTO 100
    END IF

    !! find satellite and set flag
    i=pointer_string(nprn,cprn,aprn)
    IF(i .EQ. 0) THEN
      CYCLE
    ELSE

      IF(id.EQ.1 .OR. id.EQ.2) THEN
        OB.flag(i,1)=1
        OB.lifamb(i,1,1)=0.d0
        OB.lifamb(i,1,2)=0.d0
        ambjd=0.d0
        ambjd=SIT.jdlog+(SIT.sodlog+(jepo-1)*SIT.intvlog)/86400.d0
        OB.lifamb(i,1,2)=SIT.jdlog+(SIT.sodlog+(kepo-1)*SIT.intvlog)/86400.d0
        IF (OB.lifamb(i,1,2) .GE. OB.spndel(i,1,2)) THEN
          OB.lifamb(i,1,2)=SIT.jdlog+(SIT.sodlog+(kepo-1)*SIT.intvlog)/86400.d0
          OB.lifamb(i,1,1)=MAX(ambjd,OB.spndel(i,1,2))
        ELSE
          OB.lifamb(i,1,1)=0.d0
          OB.lifamb(i,1,2)=0.d0
        END IF

        DO ifreq=2, MAXFREQ
          OB.flag(i,ifreq)=OB.flag(i,1)
          OB.lifamb(i,ifreq,1)=OB.lifamb(i,1,1)
          OB.lifamb(i,ifreq,2)=OB.lifamb(i,1,2)
        END DO

      ELSE
        deljd=0.d0
        deljd=SIT.jdlog+(SIT.sodlog+(kepo-1)*SIT.intvlog)/86400.d0
        IF (deljd .GE. OB.spndel(i,1,2)) THEN
          OB.spndel(i,1,1)=SIT.jdlog+(SIT.sodlog+(jepo-1)*SIT.intvlog)/86400.d0
          OB.spndel(i,1,2)=SIT.jdlog+(SIT.sodlog+(kepo-1)*SIT.intvlog)/86400.d0
          DO ifreq=2, MAXFREQ
            OB.spndel(i,ifreq,1)=OB.spndel(i,1,1)
            OB.spndel(i,ifreq,2)=OB.spndel(i,1,2)
          END DO
        END IF
      END IF
    END IF
  END DO

!! set obs flag
100   CONTINUE

  DO i=1,MAXSAT
    IF(OB.obs(i,1).EQ.0.d0 .OR. OB.obs(i,MAXFREQ+1).EQ.0.d0) CYCLE
    IF(((OB.jd-OB.spndel(i,1,1))*86400.d0+OB.tsec.GT.-0.5*SIT.intvlog) .AND. &
     (OB.jd-OB.spndel(i,1,2))*86400.d0+OB.tsec.LT.0.5*SIT.intvlog) THEN
      OB.obs(i,1:2*MAXFREQ)=0.d0
    END IF
  END DO

  RETURN

200 WRITE(ERROR_UNIT,'(2A)') '***ERROR(read_log_file): READ file ',LEN_TRIM(SIT.logfile)
   CALL exit(1)
  RETURN

END SUBROUTINE
