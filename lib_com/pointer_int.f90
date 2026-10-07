!*
INTEGER(IT) FUNCTION pointer_int(n,int_array,inti)
!!
!! purpose  : return the FIRST position of an element appears in an integer
!!            array.
!!            
!! parameter: n -- number of elements in the array
!!            int_array -- integer array
!!            inti      -- the element 
!!            point_int -- position of the element in the array.
!!                         zere is it is not in the array
!!
!*
USE par
IMPLICIT NONE

!*
! The arguments
!!----------------------
INTEGER(IT) n,int_array(1:*),inti

  !*
  ! The local variables
  !!-------------------------
  INTEGER(IT) i

  !*
  ! Start the exectuable code
  !!---------------------------

  pointer_int=0
  DO i=1,n
    IF (inti .EQ. int_array(i)) THEN
      pointer_int=i
      EXIT
    END IF
  END DO

  RETURN

END FUNCTION
