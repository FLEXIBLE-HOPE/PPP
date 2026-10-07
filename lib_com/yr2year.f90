!*
SUBROUTINE yr2year(iyear)
!!
!! TRANSFER THE YEAR IN 2 DIGIT TO 4 DIGIT
!!
!*
USE par
IMPLICIT NONE

!*
! THE ARGUMENTS
!!------------------------------
INTEGER(IT) :: iyear

  !*
  ! THE EXECTUABLE CODES
  !!------------------------------

  IF (iyear .GT. 1900) RETURN

  IF (iyear .LE. 79) THEN
    iyear=iyear+2000
  ELSE
    iyear=iyear+1900
  END IF

  RETURN

END SUBROUTINE


!*
SUBROUTINE year2yr(iyear)
!!
!! TRANSFER THE YEAR IN 4 DIGIT TO 2 DIGIT
!!
!*
USE par
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! THE ARGUMENTS
!!------------------------------
INTEGER(IT) :: iyear

  !*
  ! THE EXECTUABLE CODES
  !!------------------------------

  IF (iyear .GE. 2000) THEN
    iyear=iyear-2000
  ELSE IF (iyear .GE. 1945) THEN
    iyear=iyear-1900
  ELSE
    WRITE(ERROR_UNIT,'(A)') '***ERROR(year2yr): year is less 1965'
    CALL exit(1)
  END IF

  RETURN

END SUBROUTINE
