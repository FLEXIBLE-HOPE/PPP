!*
SUBROUTINE left_justify_string(string)
!!
!! TO REMOVE SPACE AT THE BEGINNING OF A STRING
!!
!! Ge Maorong (MAY.31.2003)
!!
!*
USE par
USE const
IMPLICIT NONE

!*
! THE ARGUMENTS
!!------------------------------
CHARACTER(LEN=*) :: string

  !*
  ! THE LOCAL VARIABLES
  !!------------------------------
  INTEGER(IT) :: i,j,n

  !*
  ! THE EXECTUABLE CODES
  !!------------------------------

  n=LEN(string)
  DO i=1,n
    IF (string(i:i) .NE. ' ') EXIT
  END DO
  IF (i .LE. n) string=string(i:n)

  ! EMPTY THE REST
  j=n-i+1
  DO i=j+1,n
    string(i:i)=' '
  END DO

  RETURN

END SUBROUTINE
