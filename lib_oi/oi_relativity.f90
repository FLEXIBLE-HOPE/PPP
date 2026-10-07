!*
SUBROUTINE oi_relativity(gm,gmt,xsat,xsun,acc0)
!!
!! Not complete according to IERS2003, IERS2010
!*
USE const
IMPLICIT NONE

!*
! The arguments
!!--------------------
REAL(RL) :: gm,gmt,xsat(1:*),xsun(1:*),acc0(1:*)

  !*
  ! The local variable
  !!--------------------------
  INTEGER(IT) :: i
  REAL(RL) :: r,v,rv,c1,c2,c3
  REAL(RL) :: alpha,beta,acc(3)

  REAL(RL), PARAMETER :: J(3)=[0.d0,0.d0,9.8d2]
  REAL(RL) :: crv(3),crj,cvj(3),sr(3),sv(3),csvr(3)

  !*
  ! The function called
  !!--------------------------
  REAL(RL) :: dot

  !*
  ! Start the exectable code
  !!--------------------------

  r=DSQRT(xsat(1)**2+xsat(2)**2+xsat(3)**2)
  v=DSQRT(xsat(4)**2+xsat(5)**2+xsat(6)**2)
  rv=xsat(1)*xsat(4)+xsat(2)*xsat(5)+xsat(3)*xsat(6)
  c1=gm/r
  c2=c1/r/r/(VEL_LIGHT*1.d-3)/(VEL_LIGHT*1.d-3)
  alpha=(4.d0*c1-v*v)*c2
  beta=4.d0*rv*c2

  ! Schwarzschild
  DO i=1, 3
    acc(i)=alpha*xsat(i)+beta*xsat(3+i)
  END DO

  ! Lense-Thirring precession
  CALL cross(xsat(1),xsat(4),crv)
  crj=dot(3,xsat(1),J)
  CALL cross(xsat(4),J,cvj)

  DO i=1, 3
    acc(i)=acc(i)+2.d0*c2*(3.d0*crj/r/r*crv(i)+cvj(i))
  END DO

  ! de Sitter precession
  !DO i=1, 3
  !  sr(i)=xsat(i)-xsun(i)
  !  sv(i)=xsat(i+3)-xsun(i+3)
  !END DO
  !crj=DSQRT(sr(1)**2+sr(2)**2+sr(3)**2)
  !c3=gmt/(VEL_LIGHT*1.d-3)**2/crj**3
  !DO i=1, 3
  !  sr(i)=-c3*sr(i)
  !END DO
  !CALL cross(sv,sr,csvr)
  !CALL cross(csvr,xsat(4),sv)

  ! Please note, the de Sitter has signaficant impact on POD, and the reason is
  ! unknown, as some publication, we use the scale factor 1.d-5
  !DO i=1, 3
  !  acc(i)=acc(i)+3.d0*sv(i)*1.d-5
  !END DO

  DO i=1, 3
    acc0(i)=acc0(i)+acc(i)
  END DO

  RETURN

END SUBROUTINE
