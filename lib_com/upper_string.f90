!*
CHARACTER(LEN=*) FUNCTION upper_string(string)
!!
!! LOWERCASE TO UPPERCASE
!!
!! MAORONG GE (MAY.31.2003)
!!
!*
USE par
IMPLICIT NONE

!*
! THE ARGUMENTS
!!-----------------------------
CHARACTER(LEN=*) string

  !*
  ! THE LOCAL VARIABLES
  !!-----------------------------
  INTEGER(IT) :: i,length,iadd

  !*
  ! THE EXECTUABLE CODES
  !!-----------------------------

  ! ASCII CODE DIFFERENCE BETWEEN UPPER CASE AND LOWER CASE
  iadd=ICHAR('A')-ICHAR('a')

  length=MIN(LEN(string),LEN(upper_string))
  DO i=1,length
    IF (LGE(string(i:i),'a') .AND. LLE(string(i:i),'z')) THEN
      upper_string(i:i)=CHAR(ICHAR(string(i:i))+iadd)
    ELSE
      upper_string(i:i)=string(i:i)
    END IF
  END DO

  ! LEFT THE SPACE TO THE REST
  DO i=length+1,LEN(upper_string)
    upper_string(i:i)=' '
  END DO

  RETURN

END FUNCTION
