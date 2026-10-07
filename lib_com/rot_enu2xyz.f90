!*
SUBROUTINE rot_enu2xyz(lat,lon,rotmat)
!!
!*
USE par
USE const
IMPLICIT NONE

!*
! The arguments
!!-------------------------
REAL(RL) :: lat,lon,rotmat(3,3)

  !*
  ! The local variables
  !!------------------------
  REAL(RL) :: coslat,sinlat,coslon,sinlon

  !*
  ! Start the exectuable code
  !!---------------------------

  coslat=dcos(lat-PI/2.d0)
  sinlat=dsin(lat-PI/2.d0)

  coslon=dcos(-PI/2.d0-lon)
  sinlon=dsin(-PI/2.d0-lon)

  rotmat(1,1)=coslon
  rotmat(1,2)= sinlon*coslat
  rotmat(1,3)= sinlon*sinlat
  rotmat(2,1)=-sinlon
  rotmat(2,2)=coslon*coslat
  rotmat(2,3)=coslon*sinlat
  rotmat(3,1)=0.d0
  rotmat(3,2)=-sinlat
  rotmat(3,3)=coslat
     
  RETURN

END SUBROUTINE
