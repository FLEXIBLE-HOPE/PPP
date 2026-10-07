!*
SUBROUTINE cross(v1,v2,vout)
!!
!! purpose  :  cross product of two 3-dimmension vectors
!!           
!! parameter:  v1,v2 -- 2 input 3-d vector
!!             vout  -- cross product of v1 and v2
!!
!*
USE par
IMPLICIT NONE

!*
! The arguments
!!------------------------------
REAL(RL) v1(1:*),v2(1:*),vout(1:*)

  !*
  ! The local variables
  !!---------------------
  REAL(RL) tmp(3)

  !*
  ! Start the exectuable code
  !!-----------------------------

  tmp(1)=v1(2)*v2(3)-v1(3)*v2(2)
  tmp(2)=v1(3)*v2(1)-v1(1)*v2(3)
  tmp(3)=v1(1)*v2(2)-v1(2)*v2(1)

  vout(1)=tmp(1)
  vout(2)=tmp(2)
  vout(3)=tmp(3)
 
  RETURN

END SUBROUTINE
