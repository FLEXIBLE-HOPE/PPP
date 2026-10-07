!*
SUBROUTINE read_kin_pos(flnkin,mjd,sod,snam,xcor,ikin)
!!
!*
USE par
USE const
IMPLICIT NONE

!*
! The arguments
!!--------------------------
CHARACTER(LEN=*) :: flnkin,snam
INTEGER(IT) :: mjd,ikin
REAL(RL) :: sod,xcor(3)

  !*
  ! The local variables
  !!------------------------------
  INTEGER(IT) :: lfn,mjds
  INTEGER(IT) :: iy,imon,id,ih,im

  REAL(RL) :: sec,sods,cor(3)

  CHARACTER(LEN_SITENAME) :: sitn

  LOGICAL(LG) :: lexist,lfirst
  DATA lfirst,lexist/.TRUE.,.FALSE./

  SAVE lexist,lfirst,lfn

  !*
  ! The function called
  !!------------------------------
  REAL(RL) :: timdif
  INTEGER(IT) :: modified_julday
  INTEGER(IT) :: get_valid_unit

  !*
  ! Start the exectuable code
  !!------------------------------

  IF (lfirst) THEN
    lfirst =.FALSE.
    INQUIRE(FILE=flnkin,EXIST =lexist)
    IF (lexist .EQ. .FALSE.) THEN
      ikin=-1
      RETURN
    END IF

    lfn=get_valid_unit(10)
    OPEN(FILE=flnkin,UNIT=lfn)
  END IF

  IF (lexist .EQ. .FALSE.) THEN
    ikin=-1
    RETURN
  END IF

  !! form the head, too low, need modification
  !REWIND(lfn)

  DO WHILE(.TRUE.)
    READ(lfn,'(I4,4I3,F5.1,2X,A4,3F13.3)',END=100,ERR=100) iy,imon,id,ih,im,sec,sitn,cor(1:3)
    mjds=modified_julday(id,imon,iy)
    sods=ih*3600.d0+im*60.d0+sec

    IF (DABS(timdif(mjd,sod,mjds,sods)) .LE. MAXWND) THEN
      IF (snam(1:4) .EQ. sitn) THEN
        xcor(1:3)=cor(1:3)
        ikin=0
        RETURN
      END IF
    ELSE IF (timdif(mjd,sod,mjds,sods) .LT. -MAXWND) THEN
      BACKSPACE(lfn)
      ikin=-1
      RETURN
    END IF
  END DO

100 ikin=-1
  RETURN

END SUBROUTINE
