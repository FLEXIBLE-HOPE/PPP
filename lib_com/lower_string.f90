!*
CHARACTER(LEN=*) FUNCTION lower_string(string)
!!
!! purpose  : return the lower case of the input string
!!
!! parameter: string -- input string
!!            lower_string -- return string
!!
!*
USE par
IMPLICIT NONE

!*
! The arguments
!!-----------------------
CHARACTER(LEN=*) string

  !*
  ! The local variables
  !!------------------------
  INTEGER(IT) i,length,iadd

  !*
  ! Start the exectuable code
  !!-------------------------

  !! ascii code difference between upper case and lower case
  iadd=ICHAR('a')-ICHAR('A')

  length = MIN(LEN(string),LEN(lower_string))
  DO i=1,length
    IF (LGE(string(i:i),'A') .AND. LLE(string(i:i),'Z')) THEN
      lower_string(i:i)=CHAR(ICHAR(string(i:i))+iadd)
    ELSE
      lower_string(i:i)=string(i:i)
    END IF
  END DO

  !! put SPACE for the rest
  DO i=length+1,LEN(lower_string)
    lower_string(i:i)=' '
  END DO

  RETURN

END FUNCTION
