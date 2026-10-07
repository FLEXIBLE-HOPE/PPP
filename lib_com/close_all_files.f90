!*
SUBROUTINE close_all_files()
!!
!*
USE par
IMPLICIT NONE


  !*
  ! The local variables
  !!------------------------
  INTEGER(IT) :: iunit,ierr
  CHARACTER(LEN_FILENAME):: string

  !*
  ! Start of executable code
  !!------------------------

  DO iunit=1, 300
    INQUIRE(UNIT=iunit,NAME=string)
    IF(LEN_TRIM(string) .NE. 0) THEN
      CLOSE(iunit, IOSTAT=ierr)
    END IF
  END DO

  RETURN

END SUBROUTINE
