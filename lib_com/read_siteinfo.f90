!*
SUBROUTINE read_siteinfo(SIT,jd,sod,seslen,flag)
!!
!! TO READ THE POSITION AND STOCASTIC INFORMATION FOR STATION
!!
!! GENG JIANGHUI (SEP.07.2007): CREATED
!!
!*
USE par
USE const
USE tables
USE station
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! THE ARGUMENTS
!!-----------------------------
INTEGER(IT) :: jd,flag
REAL(RL) :: sod,seslen
TYPE(SITE) :: SIT

  !*
  ! THE LOCAL VARIABLES
  !!-----------------------------
  LOGICAL(LG) :: lfirst
  INTEGER(IT) :: i,k,lfnsit,iyear,idoy,ierr
  REAL(RL) :: v(3),a(3),tref,dref,dref0
  CHARACTER(LEN_STRING) :: line, bracket, msg

  CHARACTER(LEN_FILENAME) :: sitfile

  DATA lfirst/.TRUE./
  SAVE lfirst,lfnsit

  TYPE(T_FILETABLE) :: FT

  !*
  ! THE FUNCTIONS CALLED
  !!-----------------------------
  INTEGER(IT) :: get_valid_unit
  CHARACTER(LEN_STRING) :: upper_string

  !*
  ! THE EXECTUABLE CODES
  !!-----------------------------

  IF (lfirst .EQ. .TRUE.) THEN
    lfirst=.FALSE.
    sitfile=f_tableFileName('stacoo')
    lfnsit=get_valid_unit(10)
    OPEN(UNIT=lfnsit,FILE=sitfile,STATUS='OLD',ACTION='READ',IOSTAT=ierr)
    IF (ierr .NE. 0) lfnsit=0
  END IF

  flag=-1
  SIT.x=0.d0
  SIT.dx0=0.d0
  SIT.qx=0.d0
  SIT.anttyp=''
  SIT.rectyp=''

  IF (jd .NE. 0) THEN
    CALL mjd2doy(jd,iyear,idoy)
  END IF

  IF (lfnsit .NE. 0) THEN
    REWIND(lfnsit)
    bracket='GNSS station'

    dref0=100000.d0

500 msg=upper_string(SIT.name)
    msg=msg(1:4)//'_ENU'
    DO WHILE (.TRUE.)
      READ(lfnsit,'(A)',END=200) line
      IF(INDEX(line,msg(1:8)) .NE. 0) EXIT
    END DO

    READ (line(106:),'(f18.6)') tref
    dref=(tref-(iyear*1.d0+idoy/365.25))

    IF (dref .GT. 1.d-4 ) THEN
      IF (dref .LE. dref0 ) THEN
        dref0=dref
      ELSE
        GOTO 500
      END IF
    ELSE
      IF (dref.GE.dref0 .OR. dref0.GE.0) THEN
        dref0=dref
      ELSE
        GOTO 500
      END IF
    END IF

    ! READ COORDINATES
    BACKSPACE(lfnsit)
    BACKSPACE(lfnsit)
    BACKSPACE(lfnsit)

    ! POS
    READ(lfnsit,'(A)') line
    msg=upper_string(SIT.name)
    msg=msg(1:4)//'_POS'
    IF(INDEX(line,TRIM(msg)) .NE. 0) THEN
      READ (line(10:),*,ERR=200) (SIT.x(i),i=1,3),(v(i),i=1,3),(a(i),i=1,3),(SIT.denu0(i),i=1,3)
    ELSE
      WRITE(ERROR_UNIT,'(A)') '***ERROR(read_siteinfo): station file format eror (POS)'
      WRITE(ERROR_UNIT,'(A)') TRIM(line)
    END IF

    ! SIG
    READ(lfnsit,'(A)') line
    msg=upper_string(SIT.name)
    msg=msg(1:4)//'_SIG'
    IF (INDEX(line,TRIM(msg)) .NE. 0) THEN
      READ(line(10:),'(3F14.4,3F9.4,3F9.4,3F9.4)',ERR=200) (SIT.dx0(i),i=1,9),(SIT.qx(k),k=1,3)
    ELSE
      WRITE(ERROR_UNIT,'(A)') '***ERROR(read_siteinfo): station file format eror (SIG)'
      WRITE(ERROR_UNIT,'(A)') TRIM(line)
    END IF

    READ(lfnsit,'(A)') line
    msg=upper_string(SIT.name)
    msg=msg(1:4)//'_ENU'
    IF (INDEX(line,TRIM(msg)) .NE. 0) THEN
      READ(line(10:),'(3F14.4,7X,A20,7X,A20)',ERR=200) (SIT.enu0(i),i=1,3),SIT.rectyp,SIT.anttyp
      IF (SIT.anttyp(17:20) .EQ. '    ') SIT.anttyp(17:20)='NONE'
    ELSE
      WRITE(ERROR_UNIT,'(A)') '***ERROR(read_siteinfo): station file format eror (ENU)'
      WRITE(ERROR_UNIT,'(A)') TRIM(line)
    END IF

    IF (line(132:132) .EQ. 'B') SIT.lfcb=.TRUE.

    DO i=1, 3
      SIT.x(i)=SIT.x(i)+v(i)*(iyear*1.d0+idoy/365.25-tref)
    END DO

    GOTO 500

  END IF

200 CONTINUE
  IF (SIT.x(1).NE.0.d0 .OR. SIT.x(2).NE.0.d0 .OR. SIT.x(3).NE.0.d0) THEN
    flag=0
    RETURN
  END IF

  ! THE POSITION IS NOT FOUND
  WRITE(OUTPUT_UNIT,'(A)') '%%%MESSAGE(read_sitinfo): no information for '//SIT.name

  DO i=1,3
    SIT.x(i)=1.d0
  END DO

  RETURN

ENTRY sitinfo_reset()

  lfirst=.TRUE.
  CLOSE(lfnsit)

  RETURN

END SUBROUTINE
