!*
SUBROUTINE unit_vector(n,v,u,length)
!!
!! purpose  : normalize a vector 
!!
!! parameter: n -- dimmension of  the vector
!!            v -- input vector 
!!            u -- ouput unit vector
!!            length -- length of the vector
!!           
!! author   : Ge Maorong
!!
!! last mof.: 31-May-2003 by Maorong GE, CLEAN
!*
USE par
IMPLICIT NONE

!*
! The arguments
!!---------------------
INTEGER(IT) :: n
REAL(RL) :: v(1:*),u(1:*),length

  !*
  ! The local variables
  !!--------------------
  INTEGER(IT) :: i

  length=0.d0
  DO i=1,n
    length=length+v(i)*v(i)
  END DO
  length=dsqrt(length)

  DO i=1,n
    u(i)=v(i)/length
  END DO

  RETURN

END SUBROUTINE
