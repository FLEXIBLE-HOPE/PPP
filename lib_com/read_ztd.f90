!*
SUBROUTINE read_ztd(flnztd,mjd,sitn,sztd)
!!
!*
USE par
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!---------------------
CHARACTER(LEN=*) :: flnztd
CHARACTER(LEN_SITENAME) :: sitn
REAL(RL) :: mjd,sztd

  !*
  ! The local variables
  !!--------------------------
  LOGICAL(LG) :: lfirst,lexist

  INTEGER(IT) :: i,lfn,nsit,ierr
  REAL(RL) :: ztd(MAXSIT)
  REAL(RL) :: mjd0,mjd1,t0,t1

  CHARACTER(LEN_SITENAME) :: snam(MAXSIT)
  CHARACTER(LEN_STRING) :: line

  DATA lfirst,lexist/.TRUE.,.FALSE./
  SAVE lfirst,lexist,lfn,nsit,snam,ztd,mjd0,mjd1

  !*
  ! The function called
  !!---------------------------
  INTEGER(IT) :: get_valid_unit
  INTEGER(IT) :: pointer_string

  !*
  ! Start the exectuable code
  !!----------------------------

  sztd=0.d0

  IF (lfirst .EQ. .TRUE.) THEN
    lfirst=.FALSE.

    INQUIRE(FILE=flnztd,EXIST=lexist)
    IF (lexist .EQ. .FALSE.) THEN
      WRITE(OUTPUT_UNIT,'(A)') '%%%MESSAGE(read_ztd): '//TRIM(flnztd)//' is not exist.'
      GOTO 10
    END IF

    lfn=get_valid_unit(10)
    OPEN(UNIT=lfn,FILE=flnztd)

    nsit=0
    snam=''
    ztd=0.d0
    mjd0=0.d0
    mjd1=0.d0

    line(1:1)=' '
    DO WHILE(INDEX(line,'END OF HEADER') .EQ. 0)
      READ(lfn,'(A)') line
    END DO
  END IF

10 IF (lexist .EQ. .FALSE.) THEN
    RETURN
  END IF

  IF(mjd0.EQ.0.d0 .OR. mjd.GT.mjd1) THEN

    READ(lfn,'(A)',END=200) line
    IF(line(1:3) .EQ. 'EOF') GOTO 200

    READ(line,*) snam(1),ztd(1),mjd0,mjd1

    i=1
    DO WHILE(.TRUE.)
      READ(lfn,'(A)',END=200) line
      IF(line(1:3) .EQ. 'EOF') EXIT
      IF(line(1:1) .NE. ' ') CYCLE
      i=i+1
      READ(line,*) snam(i),ztd(i),t0,t1
      IF (t0 .GT. mjd0) THEN
        i=i-1
        BACKSPACE(lfn)
        EXIT
      END IF
    END DO
    nsit=i
  END IF

  IF(mjd .GT. mjd1) GOTO 10

  i=pointer_string(nsit,snam,sitn)
  IF (i .NE. 0) THEN
    sztd=ztd(i)
  END IF

  RETURN

200 CONTINUE
  RETURN

ENTRY ztd_reset()

  mjd0=0.d0
  mjd1=0.d0
  ztd=0.d0
  nsit=0
  snam=''
  lfirst=.TRUE.
  lexist=.FALSE.
  CLOSE(lfn)

  RETURN

END SUBROUTINE
