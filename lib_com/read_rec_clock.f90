!*
SUBROUTINE read_rec_clock(flnrck,mjd,sod,sitn,rclk)
!!
!*
USE par
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!--------------------
CHARACTER(LEN=*) :: flnrck,sitn
INTEGER(IT) :: mjd
REAL(RL) :: sod, rclk

  !*
  ! The local variables
  !!---------------------------
  LOGICAL(LG) :: lfirst,lexist
  DATA lfirst,lexist/.TRUE.,.FALSE./

  INTEGER(IT) :: i,lfn,mjd0
  INTEGER(IT) :: nsit,ierr,mjdx
  INTEGER(IT) :: iy,im,id,ih,imin

  REAL(RL) :: coef(5),sec,sodx
  REAL(RL) :: a0(MAXSIT),sod0,dt

  CHARACTER(LEN_SITENAME) :: snam(MAXSIT),name
  CHARACTER(LEN_STRING) :: line

  SAVE lfn,nsit,snam,mjd0,sod0,a0,lfirst,lexist

  !*
  ! The function called
  !!-------------------------
  INTEGER(IT) :: get_valid_unit
  INTEGER(IT) :: pointer_string
  INTEGER(IT) :: modified_julday

  !*
  ! Start the exectuable code
  !!-----------------------------

  rclk=0.d0

  IF (lfirst .EQ. .TRUE.) THEN

    lfirst=.FALSE.

    INQUIRE(FILE=flnrck,EXIST=lexist)
    IF (lexist .EQ. .FALSE.) RETURN

    lfn=get_valid_unit(10)
    OPEN(UNIT=lfn,FILE=flnrck)

    i=1
    mjd0=0
    sod0=0.d0
    DO WHILE(INDEX(line,'END OF HEADER').eq.0)
      READ(lfn,'(A)',END=200) line
    END DO

  END IF

  IF(lexist .EQ. .FALSE.) RETURN

  dt=(mjd-mjd0)*86400.d0+sod-sod0

  IF(dt .LT. -MAXWND) THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(read_rec_clock): requested time before the reference time'
    CALL exit(1)

  ELSE IF(dt .GT. MAXWND) THEN

    nsit=0
    DO WHILE(.TRUE.)
      READ(lfn,'(A)',END=200,ERR=100) line
      IF(line(1:3) .NE. 'AR ') CYCLE
      coef=0.d0
      READ(line(4:),*,IOSTAT=ierr) name,iy,im,id,ih,imin,sec,i,coef(1:2)
      IF(i .EQ. 4) THEN
        READ(lfn,'(A)',END=200,ERR=100) line
        READ(line,*,IOSTAT=ierr) coef(2:3)
      ELSE IF(i .EQ. 6) THEN
        READ(lfn,'(A)',END=200,ERR=100) line
        READ(line,*,IOSTAT=ierr) coef(2:5)
      END IF

      CALL yr2year(iy)
      mjdx=modified_julday(id,im,iy)
      sodx=ih*3600.d0+imin*60.d0+sec

      dt=(mjd-mjdx)*86400.d0+(sod-sodx)
      IF(dt .LT. -MAXWND) THEN
        BACKSPACE(lfn)
        IF(i .GT. 2) BACKSPACE(lfn)
        EXIT
      ELSE IF(dt .LE. MAXWND) THEN
        IF(nsit .EQ. 0) THEN
          mjd0=mjdx
          sod0=sodx
        END IF
        nsit=nsit+1
        snam(nsit)=name
        a0(nsit)=coef(1)
      END IF
    END DO
  END IF

  i=pointer_string(nsit,snam,sitn)
  IF(i .NE. 0) THEN
    rclk=a0(i)
  END IF

  RETURN

100 CONTINUE
  WRITE(ERROR_UNIT,'(A)') '***ERROR(read_rec_clock): read '//TRIM(line)
  CALL exit(1)

200   CONTINUE
  BACKSPACE(lfn)
  IF(i .GT. 2) BACKSPACE(lfn)

  RETURN

ENTRY rec_clock_reset()

  CLOSE(lfn)
  nsit=0
  snam=''
  mjd0=0
  sod0=0.d0
  lfirst=.TRUE.
  lexist=.FALSE.

  RETURN

END SUBROUTINE
