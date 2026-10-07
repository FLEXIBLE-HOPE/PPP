!*
SUBROUTINE oi_etpole_tide(mjd,xpole,ypole,dc,ds)
!!
!*
USE par
IMPLICIT NONE

!*
! The arguments
!!--------------------
REAL(RL) :: mjd,dc,ds
REAL(RL) :: xpole,ypole

  !*
  ! The local variables
  !!---------------------------
  INTEGER(IT) :: i,j
  REAL(RL) :: xpm,ypm

  !*
  ! Start the exectuable code
  !!---------------------------

  !! m1 and m2 are in seconds of arc

  CALL mean_pole('IERS2010',mjd,xpm,ypm)
  xpm=xpole-xpm
  ypm=-(ypole-ypm)

  dc=dc-1.333d-9*(xpm+0.0115d0*ypm)
  ds=ds-1.333d-9*(ypm-0.0115d0*xpm)

  RETURN

END SUBROUTINE
