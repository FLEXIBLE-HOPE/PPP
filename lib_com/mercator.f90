!*
SUBROUTINE mercator(bl,meri,am,bm,xy)
!!
!! purpose  : transverse mercator projection
!! parameter:
!!    input : bl -- geodetic coordinates (arc)
!!            meri -- central meridian (arc)
!!            am -- semi-major axes (m)
!!            bm -- semi-minor axes (m) or inversed flatten rate
!!    output: xy -- plane coordinates (m)
!! author   : Geng J
!! created  : Nov 24, 2008
!!
!*
USE const
IMPLICIT NONE

!*
! The arguments
!!-----------------
REAL(RL) :: bl(1:*),xy(1:*),meri,am,bm

  !*
  ! The local variables
  !!-----------------------
  REAL(RL) :: a,b,e2,cme,m,m0,m2,m4,m6,m8,a0,a2,a4,a6,a8,n,lat,lon,ti,ita2,bx,arcsec

  !*
  ! Start the exectuable code
  !!-------------------------

  arcsec=180.d0*3600.d0/PI

  !! default ellipsoid WGS84
  IF (am .EQ. 0.d0) THEN
    a=6378137.d0
    b=298.257223563d0
  ELSE
    a=am
    b=bm
  END IF
  lat=bl(1)
  lon=bl(2)*180.d0/PI
  cme=meri*180.d0/PI

  !! semi axis or finv
  IF (b .LE. 6000000.d0) b=a-a/b
  e2=(a*a-b*b)/(a*a)

  n=a/DSQRT(1.d0-e2*DSIN(lat)**2)
  ti=dtan(lat)
  ita2=e2/(1.d0-e2)*DCOS(lat)**2
  m=DCOS(lat)*(lon-cme)*3600.d0/arcsec

  m0=a*(1.d0-e2)
  m2=1.5d0*e2*m0
  m4=1.25d0*e2*m2
  m6=7.d0/6.d0*e2*m4
  m8=9.d0/8.d0*e2*m6
  a0=m0+m2/2.d0+3.d0*m4/8.d0+5.d0*m6/16.d0+35.d0*m8/128.d0
  a2=m2/2.d0+m4/2.d0+15.d0*m6/32.d0+7.d0*m8/16.d0
  a4=m4/8.d0+3.d0*m6/16.d0+7.d0*m8/32.d0
  a6=m6/32.d0+m8/16.d0
  a8=m8/128.d0
  bx=a0*lat-DSIN(lat)*DCOS(lat)*(a2-a4+a6+(2.d0*a4-16.d0/3.d0*a6)*DSIN(lat)**2+16.d0/3.d0*a6*DCOS(lat)**4)

  !! result
  xy(1)=bx+n*ti*((0.5d0+(1.d0/24.d0*(5.d0-ti**2+9.d0*ita2+4.d0*ita2**2)+1.d0/72.d1*(61.d0-58.d0*ti**2+ti**4)*m**2)*m**2)*m**2)
  xy(2)=n*((1.d0+(1.d0/6.d0*(1.d0-ti**2+ita2)+1.d0/120.d0*(5.d0-18.d0*ti**2+ti**4+14.d0*ita2-58.d0*ita2*ti**2)*m**2)*m**2)*m)

  RETURN

END SUBROUTINE
