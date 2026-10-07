!*
SUBROUTINE pv2oe(xsat,oe)
!!
!! Reference:
!!           William M. Kuala. Theory of Satellite Geodesy. Blaisdell
!!           Publishing Company. 1966. P13. Equ. 2.6.
!!
!! OUTPUT:
!!         oe(1): semi-major axes (a)
!!         oe(2): eccentricities (e)
!!         oe(3): inclination (i)
!!         oe(4): Eccentric anomaly (E)
!!         oe(5): arguments of perigee (omega)
!!         oe(6): longitudes of ascending nodes (OMEGA)
!!         oe(7): argument of latitude, i.e., the sum of the true anomaly f and the angle of perigee
!*
USE const
IMPLICIT NONE

!*
! The arguments
!!---------------------
REAL(RL) :: xsat(1:*)
REAL(RL) :: oe(1:*)

  !*
  ! The local variables
  !!--------------------------
  REAL(RL) :: r,v,h,p,f,rot(3,3)
  REAL(RL) :: hr(3),pr(3),cose,sine
  
  !*
  ! Start the exectuable code
  !!--------------------------
  
  r=DSQRT(xsat(1)**2+xsat(2)**2+xsat(3)**2)
  v=DSQRT(xsat(4)**2+xsat(5)**2+xsat(6)**2)
  oe(1)=1.d0/(2.d0/r-v**2/GME)
  
  CALL cross(xsat(1),xsat(4),hr)
  h=DSQRT(hr(1)**2+hr(2)**2+hr(3)**2) 

  p=h**2/GME
  oe(2)=DSQRT(1-p/oe(1))
  
  sine=DSQRT(hr(1)**2+hr(2)**2)
  cose=hr(3)
  oe(3)=DATAN2(sine/cose,1.d0)
  
  cose=(oe(1)-r)/(oe(1)*oe(2))
  sine=r*v/(oe(2)*DSQRT(GME*oe(1)))
  oe(4)=DATAN2(sine,cose)
  
  oe(6)=DATAN2(hr(1)/(-hr(2)),1.d0)
  
  ! true anomaly
  sine=DSQRT(1-oe(2)**2)*DSIN(oe(4))
  cose=DCOS(oe(4))-oe(2)
  f=DATAN2(sine/cose,1.d0)
  
  CALL rotation(3,oe(6),xsat,pr)
  CALL rotation(1,oe(3),pr,pr)
  
  oe(7)=DATAN2(pr(2)/pr(1),1.d0)
  oe(5)=oe(7)-f
  
  RETURN

END SUBROUTINE



SUBROUTINE rotation(axis,angle,a,b)
!*
!! This subroutine computes the object of a coordinate vector while
!! rotating the coordinate system over an angle "ANGLE" around the
!! coordinate axis 'AXIS'. Vector A in the original system becomes
!! vector B in the new system. In the call to this subroutine arrays
!! A and B may be the same array.
!!
!! Reference:
!!           William M. Kuala. Theory of Satellite Geodesy. Blaisdell
!!           Publishing Company. 1966. P13. Equ. 2.6.
!!
!*
USE par
IMPLICIT NONE

!*
! Declare
!!---------------
! The ratation coordinate axis (1-X; 2-Y; 3-Z)
INTEGER(IT) :: axis
! Rotation angle (right-handed) in radians
REAL(RL) :: angle
! Vector in the original system
REAL(RL) :: a(3)
! The same vector, but now in the rotated coordinate system
REAL(RL) :: b(3)

  !*
  ! The local variables
  !!-------------------------
  INTEGER(IT) :: i1,i2,i3
  REAL(RL) :: bitwo

  !*
  ! Start the eectuable code
  !!-------------------------

  i1 = mod(axis+2,3)+1
  i2 = mod(axis+3,3)+1
  i3 = mod(axis+4,3)+1
  b(i1) = a(i1)
  bitwo = a(i2)*dcos(angle)+a(i3)*dsin(angle)
  b(i3) = a(i3)*dcos(angle)-a(i2)*dsin(angle)
  b(i2) = bitwo

  !*
  ! Return
  !!----------
  RETURN
  
END SUBROUTINE
