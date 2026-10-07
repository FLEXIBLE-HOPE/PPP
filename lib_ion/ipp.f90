!*
SUBROUTINE ipp(lat,lon,elev,azim,lat_i,lon_i)
!!
!! lat : Geodetic latitude of receiver          (rad)
!! lon : Geodetic longitude of receiver         (rad)
!! elev: Elevation angle of satellite           (rad)
!! azim: Geodetic azimuth of satellite          (rad)
!!
!*
USE const
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!----------------------
REAL(RL) :: sod,lat,lon,elev,azim,lat_i,lon_i

  !*
  ! The variables
  !!-----------------------
  REAL(RL) :: thita
  
  !*
  ! Start the exectuable code
  !!----------------------------
  
  thita = PI/2-elev-DASIN(DCOS(elev)*6.371d6/(6.371d6+4.50d5))
  
  ! lat_i in [0, pi]
  lat_i= DASIN(DSIN(lat)*DCOS(thita)+DCOS(lat)*DSIN(thita)*DCOS(azim))
  
  lon_i = lon + DASIN(DSIN(thita)*DSIN(azim)/DCOS(lat_i))

  RETURN

END SUBROUTINE