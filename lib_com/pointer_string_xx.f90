!*
INTEGER(IT) FUNCTION pointer_string_xx(n,string_array,string)
!!
!! purpose  : return the FIRST position of an element appears in a string
!!            array.
!!            
!! parameter: n -- number of elements in the array
!!            string_array --  string array
!!            string    -- the element 
!!            point_int -- position of the element in the array.
!!                         zere is it is not in the array
!!
!*
USE par
IMPLICIT NONE

!*
! The arguments
!!----------------------
INTEGER(IT) n
CHARACTER(LEN=*) string_array(1:*),string

  !*
  ! The local variables
  !!-------------------------
  INTEGER(IT) i,l

  !*
  ! Start the exectuable code
  !!--------------------------

  pointer_string_xx=0

  !! if input string is empty
  l=LEN_TRIM(string)
  IF (l .EQ. 0) RETURN

  DO i=1,n
    IF (TRIM(string_array(i)) .EQ. TRIM(string(1:l))) THEN
      pointer_string_xx=i
      EXIT
    END IF
  END DO

  RETURN

END FUNCTION
