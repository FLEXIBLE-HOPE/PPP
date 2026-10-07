!
!! purpose   : compute all tide related station position deformation (IERS2003)
!! parameters: 
!!        iers -- which terms are included
!!        jd, t  -- epoch time
!!        xsite, xsun xlun -- station-, solar- and lunar-positions
!!        rot_f2j, rot_l2f -- rotation matrix from earth-fixed to inertial
!!                and from station to earth-fixed system
!!        xpole,ypole -- sideral time, x and y pole positions
!!        olc  -- ocean loading coefficients
!!        dx  -- position correction
!
SUBROUTINE tide_displace(mjd,sod,mjdut1,sodut1,xsit_f,xsun,xlun,rot_f2j,rot_l2f,&
                         lat,lon,xpole,ypole,olc,disp)
USE par
USE const
IMPLICIT NONE

INTEGER(IT) :: mjd,mjdut1
REAL(RL) :: sod,sodut1,xsit_f(1:*),xsun(1:*),xlun(1:*),disp(1:*)
REAL(RL) :: lat,lon,xpole,ypole,rot_f2j(3,3),rot_l2f(3,3),olc(11,6)

  !*
  ! The local variables
  !!---------------------
  INTEGER(IT) :: i,j,mjdutc,mjdtai
  REAL(RL) :: xpm,ypm,sodutc,sodtai,dxi(3),colat,xs(3),xl(3),dump(2)

  !*
  ! The function called
  !!---------------------
  REAL(RL) :: dot

  !*
  ! Start the exectuable code
  !!--------------------------

  !! initialization
  DO i=1,3
    disp(i)=0.d0
  END DO
  colat=PI/2.D0-lat

  CALL timinc(mjd,sod,OFF_GPS2TAI,mjdtai,sodtai)
  !CALL iau_TAIUTC(DBLE(mjdtai+2400000.5D0),sodtai/86400.d0,dump(1),dump(2),i)
  !mjdutc=INT(dump(1)-2400000.5d0+dump(2))
  !sodutc=(dump(1)-2400000.5d0+dump(2)-mjdutc)*86400.d0
  !IF (DABS(sodutc-NINT(sodutc)) .LE. 5.d-5) sodutc=DBLE(NINT(sodutc))
  CALL taiutc(mjdtai,sodtai,mjdutc,sodutc)

  !! 1. Displacement due to frequency-independent solid-Earth tide(in J2000)
  !! The tidal model contains a time-independent part so that the coordinates obtained by taking
  !! into account this model in the analysis will be "conventional tide free" values.
  CALL matmpy(xsun,rot_f2j,xs,1,3,3)
  CALL matmpy(xlun,rot_f2j,xl,1,3,3)
  CALL DEHANTTIDEINEL(xsit_f(1:3),mjdutc,sodutc/3600,xs,xl,dxi)
  CALL matmpy(rot_f2j,dxi,dxi,3,3,1)
  DO i=1,3
    disp(i)=disp(i)+dxi(i)
  END DO

  !! 1.1 Permanent deformation ENU
  !! compute "mean tide" coordinates from "conventional tide free" coordinates
  !! dxi(2)=-0.0252-0.0001*(3*DSIN(lat)**2-1.d0)/2.d0
  !! dxi(3)=(-0.1206+0.0001*(3*DSIN(lat)**2-1.d0))*(3*DSIN(lat)**2-1.d0)
  !! disp(1:3)=disp(1:3)+dxi(1:3)

  !! 2. Displacement due to the pole motion
  !! xpole,ypole in seconds of arc (equ.22, pp 67)
  !! dxi(3) east-north-radial
  CALL mean_pole('IERS2010',mjd+(sod+OFF_GPS2TT)/86400.d0,xpm,ypm)

  xpm=xpole-xpm
  ypm=-(ypole-ypm)
  ! east
  dxi(1)=  9.d0*dcos(colat)     *(xpm*dsin(lon)-ypm*dcos(lon))
  ! from south to north
  dxi(2)=  9.d0*dcos(2.d0*colat)*(xpm*dcos(lon)+ypm*dsin(lon))
  ! upword
  !! IERS 2003
  !!dxi(3)=-32.d0*dsin(2.d0*colat)*(xpm*dcos(lon)+ypm*dsin(lon))
  !! IERS 2010
  dxi(3)=-33.d0*dsin(2.d0*colat)*(xpm*dcos(lon)+ypm*dsin(lon))

  !! rotation matrix from east-north-radial to x-y-z, then to J2000
  CALL matmpy(rot_l2f,dxi,dxi,3,3,1)
  CALL matmpy(rot_f2j,dxi,dxi,3,3,1)
  DO i=1,3
    disp(i)=disp(i)+dxi(i)*1.d-3
  END DO

  !! 3. Displacement due to ocean-loading
  dxi(1:3)=0.d0
  !! CALL ocean_tidal_loading(mjdutc,sodutc,olc,1,0.d0,dxi)
  CALL HARDISP(mjdutc,sodutc,olc,1,0.d0,dxi)
  CALL matmpy(rot_l2f,dxi,dxi,3,3,1)
  CALL matmpy(rot_f2j,dxi,dxi,3,3,1)
  DO i=1,3
    disp(i)=disp(i)+dxi(i)
  END DO

  !! ocean pole tide loading

  !! S1-S2 atmospheric pressure loading
  !dxi(1:3)=0.d0
  !! the time should be ut1, but we use utc to replace ut1 for similaty
  !CALL GRDINTRP(mjdut1+sodut1/86400.d0,lat,lon,dxi)
  !CALL matmpy(rot_l2f,dxi,dxi,3,3,1)
  !CALL matmpy(rot_f2j,dxi,dxi,3,3,1)
  !DO i=1,3
  !  disp(i)=disp(i)+dxi(i)
  !END DO

  RETURN

END SUBROUTINE
