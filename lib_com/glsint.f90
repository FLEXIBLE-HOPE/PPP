!*
SUBROUTINE glsrkf4(h,x,eph)
!!
!*
USE brdeph
IMPLICIT NONE

!*
! The arguments
!!------------------------
REAL(RL) :: h,x(1:*)
TYPE(GLONASS_BRDEPH) :: eph

  !*
  ! The local variables
  !!--------------------------
  INTEGER(IT) :: i,j
  REAL(RL) :: coef(4),funct(6),acc(6),xj(6)

  !*
  ! Start the exectuable code
  !!--------------------------

  coef(1) = h/2.0
  coef(2) = h/2.0
  coef(3) = h
  coef(4) = h

  DO i=1, 6
    xj(i) = x(i)
  END DO

  CALL glsfright(xj,acc,eph)

  DO i=1, 6
    funct(i)=Xj(i)
  END DO

  DO j=1, 3
    DO i=1, 6
       xj(i)=x(i)+coef(j)*acc(i)
       funct(i)=funct(i)+coef(j+1)*acc(i)/3.0d0
    END DO
    CALL glsfright(xj,acc,eph)
  END DO

  DO i=1, 6
    x(i)=funct(i)+h*acc(i)/6.0d0
  END DO

  RETURN

END SUBROUTINE


!*
SUBROUTINE glsfright(x,acc,eph)
!!
!*
USE brdeph
IMPLICIT NONE

!*
! The arguments
!!------------------------
REAL(RL) :: x(1:*),acc(1:*)
TYPE(GLONASS_BRDEPH) :: eph

  !*
  ! Local variables
  !!-----------------
  INTEGER(IT) :: i
  REAL(RL) :: vec1(3), vec2(3), r
  REAL(RL), PARAMETER :: C20 = -1.08262575D-3
  REAL(RL), PARAMETER :: GMEARTH = 3.986004418D5
  REAL(RL), PARAMETER :: EROT = 7.292115D-5
  REAL(RL), PARAMETER :: ER = 6.378136D3
  REAL(RL) :: factor1, factor2

  !*
  ! Start the exectuable code
  !!--------------------------------

  r = DSQRT(x(1)**2+x(2)**2+x(3)**2)
  factor1 = -GMEARTH/(r**3)
  factor2 = 1.5d0*GMEARTH*C20*ER*ER/(r**5)

  vec1(1) = EROT*EROT*X(1)+2.0d0*EROT*X(5)
  vec1(2) = EROT*EROT*X(2)-2.0d0*EROT*X(4)
  vec1(3) = 0.d0

  vec2(1) = x(1)*(1.0d0-5.0d0*((x(3)/r)**2))
  vec2(2) = x(2)*(1.0d0-5.0d0*((x(3)/r)**2))
  vec2(3) = x(3)*(3.0d0-5.0d0*((x(3)/r)**2))

  DO i=1, 3
    acc(i) = X(i+3)
    acc(i+3) = eph.acc(i)+vec1(i)+factor2*vec2(i)+factor1*X(i)
  END DO

  RETURN

ENTRY glsinit(x,eph)

  DO i=1, 3
    x(i)=eph.pos(i)
  END DO

  DO i=4, 6
    x(i)=eph.vel(i-3)
  END DO

  RETURN

END SUBROUTINE
