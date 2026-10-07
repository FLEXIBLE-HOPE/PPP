!
SUBROUTINE oi_aod(mjd,lmax,dca,dsa)
!!
!*
USE tables
IMPLICIT NONE

!*
! The arguments
!!--------------------
REAL(RL) :: mjd
INTEGER(IT) :: lmax
REAL(RL) :: dca(MAXAODDEG,0:MAXAODDEG)
REAL(RL) :: dsa(MAXAODDEG,0:MAXAODDEG)


  !*
  ! The local variables
  !!------------------------------
  TYPE AOD_MOD
    REAL(RL) :: mjd
    REAL(RL) :: c(0:MAXAODDEG,0:MAXAODDEG)
    REAL(RL) :: s(0:MAXAODDEG,0:MAXAODDEG)
  END TYPE
  TYPE(AOD_MOD) :: r(2)

  LOGICAL(LG) :: lfirst,lexist,lfound
  INTEGER(IT) :: i,j,lfn,l,m,l_in,m_in
  INTEGER(IT) :: iy,imon,id,ih,im
  REAL(RL) :: mjd0,mjds,mjde,sod,alpha,is

  CHARACTER(LEN_STRING) :: line,flnaod

  DATA lfirst/.TRUE./
  SAVE lfirst,lfn,mjds,mjde,r

  !*
  ! The function called
  !!------------------------------
  INTEGER(IT) :: get_valid_unit

  !*
  ! Start the exectuable code
  !!------------------------------

  IF (lfirst .EQ. .TRUE.) THEN
    lfirst=.FALSE.
    !! TT to GPS
    mjd0=mjd-OFF_GPS2TT/86400.d0
    CALL mjd2date(INT(mjd0),(mjd0-INT(mjd0))*86400.d0,iy,imon,id,ih,im,is)

    CALL get_file_name(.FALSE.,'aod',' ',iy,imon,id,ih,flnaod)
    INQUIRE(FILE=flnaod,EXIST=lexist)
    IF (lexist .EQ. .FALSE.) THEN
      WRITE(OUTPUT_UNIT,'(A)') '###WARNING(oi_aod): '//TRIM(flnaod)//' is not exist'
      RETURN
    END IF
    lfn=get_valid_unit(10)
    OPEN(UNIT=lfn,FILE=flnaod)

    READ(lfn,'(A/A/A)',END=100) line,line,line
    IF (INDEX(line,'%% Start and stop time :') .EQ. 0) THEN
      WRITE(OUTPUT_UNIT,'(A)') '***ERROR(oi_aod): The third line does not contain start&stop time.'
      WRITE(OUTPUT_UNIT,'(A)') TRIM(line)
      CALL exit(1)
    END IF
    READ(line(25:),*,ERR=100) i,mjds,j,mjde
    mjds=i+mjds/86400.d0
    mjde=j+mjde/86400.d0

    READ(lfn,*) j,sod
    r(1).mjd=j+sod/86400.d0
    ! Reading the first model
    DO l=0, MAXAODDEG
      DO m=0, l
        READ(lfn,*,ERR=100) l_in,m_in,r(1).c(l,m),r(1).s(l,m)
        IF (l_in.NE.l .OR. m_in.NE.m) GOTO 100
      END DO
    END DO

    READ(lfn,*) j,sod
    r(2).mjd=j+sod/86400.d0
    ! Reading the second model
    DO l=0, MAXAODDEG
      DO m=0, l
        READ(lfn,*,ERR=100) l_in,m_in,r(2).c(l,m),r(2).s(l,m)
        IF (l_in.NE.l .OR. m_in.NE.m) GOTO 100
      END DO
    END DO
  END IF

  !! check data span
  mjd0=mjd-OFF_GPS2TT/86400.d0
  IF (mjd0.LT.mjds .OR. mjd0 .GT. mjde) RETURN

  !! check the current record
  lfound=.FALSE.
  DO WHILE(lfound .EQ. .FALSE.)
    IF (mjd0.GE.r(1).mjd-MAXWND/86400.d0 .AND. mjd0.LE.r(2).mjd+MAXWND/86400.d0) THEN
      lfound=.TRUE.
    ELSE IF (mjd0 .GT. r(2).mjd) THEN
      r(1)=r(2)
      READ(lfn,*) j,sod
      r(2).mjd=j+sod/86400.d0
      DO l=0, MAXAODDEG
        DO m=0,l
          READ(lfn,*,ERR=100) l_in,m_in,r(2).c(l,m),r(2).s(l,m)
          IF (l_in.NE.l .OR. m_in.NE.m) GOTO 100
        END DO
      END DO
    ELSE
      r(2)=r(1)
      BACKSPACE(lfn)
      DO l=0, MAXAODDEG
        DO m=0,l
          BACKSPACE(lfn)
        END DO
      END DO
      BACKSPACE(lfn)
      DO l=0, MAXAODDEG
        DO m=0,l
          BACKSPACE(lfn)
        END DO
      END DO
      READ(lfn,*) j,sod
      r(1).mjd=j+sod/86400.d0
      DO l=0, MAXAODDEG
        DO m=0,l
          READ(lfn,*,ERR=100) l_in,m_in,r(1).c(l,m),r(1).s(l,m)
          IF (l_in.NE.l .OR. m_in.NE.m) GOTO 100
        END DO
      END DO
    END IF
  END DO

  !! not found
  IF (lfound .EQ. .TRUE.) THEN
    WRITE(OUTPUT_UNIT,'(A)') '***ERROR(oi_aod): requested time not found,',mjd0,r(1).mjd,r(2).mjd
    CALL exit(1)
  END IF

  !! interpolation
  alpha=(mjd0-r(1).mjd)/(r(2).mjd-r(1).mjd)

  ! degree 0 and 1 are excluded
  DO l=2, lmax
    DO m=0, l
      dca(l,m) = (1.d0-alpha)*r(1).c(l,m)+alpha*r(2).c(l,m)
      dsa(l,m) = (1.d0-alpha)*r(1).s(l,m)+alpha*r(2).s(l,m)
    END DO
  END DO

100  RETURN

END SUBROUTINE
