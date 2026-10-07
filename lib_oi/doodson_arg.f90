!*
SUBROUTINE doodson_arg(gmst,arg,beta)
!!
!*
USE const
IMPLICIT NONE

!*
! The arguments
!!-----------------------------
! Greenwich Mean Sidereal Time [radian]
REAL(RL) :: gmst
! fundamental arguments of nutation theory
REAL(RL) :: arg(1:*)
! doodsnon's fundamental arguments
REAL(RL) :: beta(1:*)

  !*
  ! Start the exectuable code
  !!---------------------------

  beta(2)=arg(3)+arg(5)
  beta(3)=beta(2)-arg(4)
  beta(4)=beta(2)-arg(1)
  beta(5)=-arg(5)
  beta(6)=beta(2)-arg(4)-arg(2)
  beta(1)=gmst+PI-beta(2)

  RETURN

END SUBROUTINE
