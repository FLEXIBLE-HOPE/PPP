!*
SUBROUTINE ion2(mjd,sod,xsit,xsat,stec)
!!
!* 
USE const
IMPLICIT NONE

!*
! The arguments
!!---------------------------
REAL(RL) :: xsit(3)
REAL(RL) :: xsat(3)
REAL(RL) :: sod,stec
INTEGER(IT) :: mjd

  !*
  ! The local variables
  !!---------------------------
  REAL(RL):: bcos,year
  INTEGER(IT) :: iy,idoy,is,ie

  INTEGER(IT) :: modified_julday

  !*
  ! Start the exectuable code
  !!---------------------------

  CALL mjd2doy(mjd,iy,idoy)
  
  is=modified_julday(1,1,iy)
  ie=modified_julday(31,12,iy)

  year=(idoy+sod/86400.d0)/(ie-is+1)+iy

  CALL pathIntegral(year,xsit,xsat,bcos)

  stec=stec*1.0d7*bcos

  RETURN

END SUBROUTINE


!*
SUBROUTINE PathIntegral(year,xsit,xsat,itgr)
!!
!!
USE const
IMPLICIT NONE

!*
! The arguments
!!---------------------------
REAL(RL) :: year,itgr
REAL(RL) :: xsit(3)
REAL(RL) :: xsat(3)

  !*
  ! The local variables
  !!---------------------------
  REAL(RL) :: dx,dy,dz,point(3),sca,neu(3)
  REAL(RL) :: geod(3),f,length
  REAL(RL) :: rot_l2f(3,3)

  !*
  ! Start the exectuable code
  !!---------------------------

  ! in meter
  dx=xsat(1)-xsit(1)
  dy=xsat(2)-xsit(2)
  dz=xsat(3)-xsit(3)

  sca=400000.d0/DSQRT(dx*dx+dy*dy+dz*dz)

  point(1)=xsit(1)+sca*dx
  point(2)=xsit(2)+sca*dy
  point(3)=xsit(3)+sca*dz

  CALL xyzblh(point,1.D0,0.D0,0.D0,0.D0,0.D0,0.D0,geod)

  CALL igrf13syn(0,year,1,400.d0,90.0-geod(1)*RAD2DEG,geod(2)*RAD2DEG,dx,dy,dz,f)

  neu(1)=xsat(1)-point(1)
  neu(2)=xsat(2)-point(2)
  neu(3)=xsat(3)-point(3)
  CALL rot_enu2xyz(geod(1),geod(2),rot_l2f)
  CALL matmpy(neu,rot_l2f,neu,1,3,3)

  length=(dx*neu(1)+dy*neu(2)-dz*neu(3))/(DSQRT(neu(1)*neu(1)+neu(2)*neu(2)+neu(3)*neu(3))*f)

  itgr=DABS(length)*f

  RETURN

END SUBROUTINE
