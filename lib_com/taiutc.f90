!*
SUBROUTINE taiutc(mjdtai,sodtai,mjdutc,sodutc)
!!
!*
USE par
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!-------------------
INTEGER(IT) :: mjdtai,mjdutc
REAL(RL) :: sodtai,sodutc

  !*
  ! The local variables
  !!----------------------
  INTEGER(IT) :: i,mjd0,mjd1
  REAL(RL) :: sod0,sod1,timdif

  !*
  ! Start the exectuable code
  !!--------------------------

  ! Iterate (though in most cases just once is enough)
  mjd0=mjdtai
  sod0=sodtai
  DO i=1, 3

    ! Guessed UTC to TAI.
    CALL utctai(mjd0,sod0,mjd1,sod1)

    ! Adjust guessed UTC.
    CALL timinc(mjd0,sod0,timdif(mjdtai,sodtai,mjd1,sod1),mjd0,sod0)
  END DO

  mjdutc=mjd0
  sodutc=sod0

  RETURN

END SUBROUTINE


SUBROUTINE utctai(mjdutc,sodutc,mjdtai,sodtai)
!!
!*
USE par
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!-------------------
INTEGER(IT) :: mjdtai,mjdutc
REAL(RL) :: sodtai,sodutc

  !*
  ! The local variables
  !!----------------------
  REAL(RL) :: ls0,ls24
  REAL(RL) :: dleap,sod

  !*
  ! The function called
  !!-----------------------
  REAL(RL) :: leap

  !*
  ! Start the exectuable code
  !!--------------------------

  ! Get TAI-UTC at 0h today
  ls0=leap(mjdutc)

  ! Get TAI-UTC at 0h tomorrow (to detect jumps).
  ls24=leap(mjdutc+1)

  dleap=ls24-ls0

  ! Remove any scaling applied to spread leap into preceding day.
  sod=sodutc*dleap/86400.d0

  CALL timinc(mjdutc,sodutc,sod+ls0,mjdtai,sodtai)

  RETURN

END SUBROUTINE


REAL(RL) FUNCTION leap(mjdutc)
!!
!*
USE par
USE tables
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!-------------------
INTEGER(IT) :: mjdutc

  !*
  ! The local variables
  !!----------------------
  INTEGER(IT) :: lun,ios,i
  LOGICAL(LG) :: first
  INTEGER(IT) :: jdt(50),ls(50),nls
  CHARACTER(LEN_STRING) :: line
  CHARACTER(LEN_FILENAME) :: leapfile
  TYPE(T_FILETABLE) :: FT

  DATA first /.TRUE./
  SAVE nls,jdt,ls,first

  !*
  ! The funciton called
  !!----------------------
  INTEGER(IT) :: get_valid_unit

  !*
  ! Start the exectuable code
  !!-----------------------------

  IF (first .EQ. .TRUE.) THEN
    leapfile = f_tablefilename('leapsc')
    lun=get_valid_unit(10)
    OPEN(UNIT=lun,FILE=leapfile,STATUS='old',IOSTAT=ios)
    IF (ios .NE. 0) THEN
      WRITE(ERROR_UNIT,'(A,A)') '***ERROR(leap): open ', TRIM(leapfile)
      CALL exit(1)
    END IF
    first=.FALSE.

    !! first line
    DO WHILE(INDEX(line,'+leap sec') .EQ. 0)
      READ(lun,'(A)',END=100) line
    END DO

    !! read file
    nls=0
    READ(lun,'(A)',END=100) line
    DO WHILE(INDEX(line,'-leap sec') .EQ. 0)
      nls=nls+1
      IF (nls .GT. 50) THEN
        WRITE(ERROR_UNIT,'(A)') '***ERROR(leap): more than 50 LS records, change dimmension'
        CALL exit(1)
      END IF
      READ(line,*,ERR=200)  jdt(nls),ls(nls)
      READ(lun,'(A)',END=100) line
    END DO
    CLOSE(lun)
  END IF

  IF (mjdutc .LE. jdt(1)) THEN
    WRITE(ERROR_UNIT,'(A,2I7)') '***ERROR(leap): epoch before table start,',mjdutc,jdt(1)
    CALL exit(1)
  END IF

  IF (mjdutc .GT. jdt(nls)) THEN
    WRITE(ERROR_UNIT,'(A,2I7)') '***ERROR(leap): epoch after table end,',mjdutc,jdt(nls)
    CALL exit(1)
  END IF

  DO i=1,nls-1
    !IF (mjdutc.GT.jdt(i) .AND. mjdutc.LE.jdt(i+1)) THEN
    IF (mjdutc.GE.jdt(i) .AND. mjdutc.LT.jdt(i+1)) THEN
      leap=ls(i)
      RETURN
    END IF
  END DO

  RETURN

100 CONTINUE
  WRITE(ERROR_UNIT,'(A)') '***ERROR(leap): sinex endline -leap sec not found'
  CALL exit(1)

200 CONTINUE
  WRITE(ERROR_UNIT,'(A/A)') '***ERROR(leap): read leap.sec error,',line
  CALL exit(1)

END FUNCTION
