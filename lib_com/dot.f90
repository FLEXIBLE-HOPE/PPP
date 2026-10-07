!*
REAL(RL) FUNCTION dot(n,v1,v2)
!!
!! purpose  : dot product of twe vectors
!!
!! parameter: n -- dimmesion of vectors
!!             v1, v2 -- input vecotrs
!!             dot -- dot product of v1 and v2
!*
USE par
IMPLICIT NONE

!*
! The arguments
!!------------------------
INTEGER(IT) n
REAL(RL) v1(1:*),v2(1:*)

  !*
  ! The local variables
  !!----------------------
  INTEGER(IT) i

  !*
  ! Start the exectuable code
  !!---------------------------

  dot = 0.d0
  DO i=1,n
    dot=dot+v1(i)*v2(i)
  END DO

  RETURN

END FUNCTION
