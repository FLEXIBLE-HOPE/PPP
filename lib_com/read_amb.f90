!*
SUBROUTINE read_amb(flnamb,mjd,snam,cprn,xamb)
!!
!*
USE par
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!-------------------
CHARACTER(LEN=*) :: flnamb
CHARACTER(LEN_SITENAME) :: snam
CHARACTER(LEN_PRN) :: cprn
REAL(RL) :: mjd,xamb

  !*
  ! The local variables
  !!----------------------------

  INTEGER(IT) :: lfn

  REAL(RL) :: t0,t1,xini,xcor

  CHARACTER(LEN_PRN) :: prn
  CHARACTER(LEN_SITENAME) :: sitn
  CHARACTER(LEN_STRING) :: line

  LOGICAL(LG) :: lfirst,lexist
  DATA lfirst /.TRUE./

  SAVE lfirst,lfn

  !*
  ! The function used
  !!----------------------------
  INTEGER(IT) :: get_valid_unit

  !*
  ! Start the exectuable code
  !!----------------------------

  IF (lfirst .EQ. .TRUE.) THEN
    lfirst=.FALSE.

    INQUIRE(FILE=flnamb,EXIST=lexist)
    IF (lexist .EQ. .FALSE.) THEN
      WRITE(OUTPUT_UNIT,'(A)') '%%%MESSAGE(read_amb): '//TRIM(flnamb)//' is not exist.'
      CALL exit(1)
    END IF

    lfn=get_valid_unit(10)
    OPEN(UNIT=lfn,FILE=flnamb)

  END IF

  xamb=0.d0

  REWIND(lfn)

  DO WHILE(.TRUE.)
    READ(lfn,'(A)',END=200) line
    IF (line(8:10) .NE. 'AMB') CYCLE

    READ(line,'(23X,A4,1X,A3,2F25.10,2F18.10)',END=200) sitn,prn,xini,xcor,t0,t1

    IF (sitn .NE. snam) CYCLE
    IF (prn .NE. cprn) CYCLE

    IF (mjd.GT.(t0-MAXWND/86400.d0) .AND. mjd.LT.(t1+MAXWND/86400.d0)) THEN
      xamb=xini+xcor
      EXIT
    END IF

  END DO

200 CONTINUE

  RETURN

END SUBROUTINE
