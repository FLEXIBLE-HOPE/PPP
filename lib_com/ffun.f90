!*
REAL(RL) FUNCTION ffun(phi,h)
!!
!! Computes the ellipsoidal/elevation-dependent variation of the
!! average acceleration of gravity from Saastamoinen model.
!!
!! INPUT:
!!  PHI    Geocentric latitude, radians
!!  H      Elevation of site above geoid, km
!*
USE par
IMPLICIT NONE

!*
! The arguments
!!-------------------
REAL(RL) :: phi,h

  !*
  ! Start the exectuable code
  !!----------------------------------
  ffun=1.d0-0.266D-2*DCOS(2.d0*phi)-0.28d-3*h

  RETURN

END FUNCTION
