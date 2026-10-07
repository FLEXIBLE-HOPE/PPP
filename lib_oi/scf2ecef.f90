!!
SUBROUTINE scf2ecef(x,r)
!!
!*
USE par
IMPLICIT NONE

!*
! The arguments
!!------------------------
REAL(RL) :: x(1:*),r(3,3)

  !*
  ! The local variables
  !!----------------------------
  INTEGER(IT) :: i
  REAL(RL) :: ux(3),uv(3),un(3),lent

  !*
  ! Start the exectuable code
  !!----------------------------

  DO i=1, 3
    ux(i)=x(i)
    uv(i)=x(3+i)
  END DO

  CALL unit_vector(3,ux,ux,lent)
  CALL unit_vector(3,uv,uv,lent)
  CALL cross(ux,uv,un)
  CALL unit_vector(3,un,un,lent)
  ! previous using this for GNSS and scf for LEO
  !CALL cross(un,ux,uv)
  !CALL unit_vector(3,uv,uv,lent)
  CALL cross(uv,un,ux)
  CALL unit_vector(3,ux,ux,lent)

  DO i=1, 3
    r(i,1)=uv(i)
    r(i,2)=un(i)
    r(i,3)=ux(i)
  END DO

  RETURN

END SUBROUTINE
