!*
SUBROUTINE xyzblh(x,scale,a0,b0,dx,dy,dz,geod)
!!
!! purpose  :  computation of ellipsoidal coordinates rb,rl,rh
!!             given the cartesian coordinates x,y,z
!!
!! parameter:    x,y,z:    cartesian coord. of the point in metres
!!               scale:    scale to meter
!!            dx,dy,dz:    translation components from the origin of
!!                         the cart. coord. system (x,y,z) to the
!!                         center of the ref.ellipsoid (in metres)
!!                 a,b:    semi-major and semi-minor axes of the
!!                         ref.ellipsoid in metres
!!                  rb:    ell. latitude (arc)
!!                  rl:    ell. longitude (pos. east of greenwich)
!!                  rh:    ell. height (m)
!!
!! modified :  default ellipsoid is WGS84
!!
!*
USE par
USE const
IMPLICIT NONE

!*
! THE ARGUMENTS
!!-----------------------------
REAL(RL) :: x(1:*),scale,a0,b0,dx,dy,dz,geod(1:*)

  !*
  ! THE LOCAL VARIABLES
  !!-----------------------------
  REAL(RL) :: a,b,xp,yp,zp,e2,s,zps,n,hp,rbi,rbd,rhd

  !*
  ! THE EXECTUABLE CODES
  !!-----------------------------

  ! THE DEFAULT ELLIPSOID IS WGS84
  IF (a0 .EQ.0.d0) THEN
    a=6378137.d0
    b=298.257223563d0
  ELSE
    a=a0
    b=b0
  END IF

  ! TO METERS
  xp=x(1)*scale+dx
  yp=x(2)*scale+dy
  zp=x(3)*scale+dz

  ! SEMI-AXIS
  IF (b .LE. 6000000.d0) b=a-a/b
  e2=(a*a-b*b)/(a*a)

  s=DSQRT(xp*xp+yp*yp)
  geod(2)=DATAN2(yp,xp)
  IF (geod(2) .LT. 0.d0) geod(2)=2.d0*PI+geod(2)
  zps=zp/s
  geod(3)=DSQRT(xp*xp+yp*yp+zp*zp)-a
  geod(1)=DATAN(zps/(1.d0-e2*a/(a+geod(3))))

10 n=a/DSQRT(1.d0-e2*DSIN(geod(1))**2)
  hp=geod(3)
  rbi=geod(1)
  geod(3)=s/DCOS(geod(1))-n
  geod(1)=DATAN(zps/(1.d0-e2*n/(n+geod(3))))
  rbd=DABS(rbi-geod(1))
  rhd=DABS(hp-geod(3))
  IF (rbd*n.GT.1.d-4 .OR. rhd.GT.1.d-4) GOTO 10

  RETURN

 END SUBROUTINE
