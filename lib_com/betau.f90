!*
SUBROUTINE betau(xsat,xsun,beta,u)
!!
!*
USE const
IMPLICIT NONE

!*
! The arguments
!!-------------------
REAL(RL) :: beta,u
REAL(RL) :: xsat(1:*),xsun(1:*)

  !*
  ! The local variables
  !!--------------------------
  REAL(RL) :: xsat_unit(3),nop_unit(3)
  REAL(RL) :: xsun_unit(3),vsat_unit(3),det

  !*
  ! The function called
  !!--------------------------
  REAL(RL) :: dot

  !*
  ! Start the exectuable code
  !!--------------------------

  CALL unit_vector(3,xsat,xsat_unit,det)
  CALL cross(xsat,xsat(4),nop_unit)
  CALL unit_vector(3,nop_unit,nop_unit,det)
  CALL unit_vector(3,xsun,xsun_unit,det)
  beta=PI/2.d0-DACOS(dot(3,nop_unit,xsun_unit))

  u=DACOS(dot(3,xsat_unit,xsun_unit)/DCOS(beta))
  CALL unit_vector(3,xsat(4),vsat_unit,det)
  IF (dot(3,vsat_unit,xsun_unit) .GT. dot(3,vsat_unit,xsat_unit)*dot(3,xsat_unit,xsun_unit)) u=-u
  u=u-PI

  DO WHILE(u .LT. 0.0)
    u=u+2*PI
  END DO
  DO WHILE(u .GT. 2*PI)
    u=u-2*PI
  END DO

  RETURN

END SUBROUTINE
