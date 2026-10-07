!*
INTEGER(IT) FUNCTION get_valid_unit(iunit0)
!!
!! FINDING A VALID FILE UNIT FOR A NEW FILE
!!
!*
USE par
USE const
IMPLICIT NONE

  !*
  ! THE ARGUMENTS
  !!------------------------------
  INTEGER(IT) :: iunit0

  !*
  ! THE LOCAL VARIABLES
  !!------------------------------
  INTEGER(IT) :: iunit
  CHARACTER(LEN_FILENAME) :: string

  !*
  ! THE EXECTUABLE CODE
  !!------------------------------

  get_valid_unit=0
  iunit=iunit0
  IF (iunit0 .LE. 10) iunit=10
  DO WHILE(get_valid_unit .EQ. 0)
    INQUIRE(UNIT=iunit,NAME=string)
    IF (LEN_TRIM(string) .NE. 0) THEN
      iunit=iunit+1
    ELSE
      get_valid_unit=iunit
    END IF
  END DO

  RETURN

END FUNCTION


!*
INTEGER(IT) FUNCTION get_file_unit(fname)
!!
!*
USE par
IMPLICIT NONE

!*
! THE ARGUMENT
!!------------------------------
CHARACTER(LEN=*) :: fname

  !*
  ! THE LOCAL VARIABLES
  !!------------------------------
  INTEGER(IT) :: lfn
  CHARACTER(LEN_FILENAME):: string

  !*
  ! THE EXECTUABLE CODE
  !!------------------------------

  get_file_unit=0

  lfn=10
  DO WHILE (get_file_unit.EQ.0 .AND. lfn.LE.1000)
    INQUIRE(UNIT=lfn,NAME=string)
    string=string(MAX(INDEX(string,'\',.TRUE.),INDEX(string,'/',.TRUE.))+1:)
    IF (string(1:LEN_TRIM(string)) .EQ. fname(1:LEN_TRIM(fname))) THEN
      get_file_unit=lfn
      RETURN
    END IF
    lfn=lfn+1
  END DO

  RETURN

END FUNCTION
