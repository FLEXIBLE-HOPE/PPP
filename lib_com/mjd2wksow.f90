!*
SUBROUTINE mjd2wksow(mjd,sod,week,sow)
!*
USE par
IMPLICIT NONE

!*
! The arguments
!!------------------------
INTEGER(IT) :: mjd, week
REAL(RL) :: sod, sow

  !*
  ! Start the exectuable code
  !!----------------------------

  week=INT((mjd+sod/86400.d0-44244d0)/7)
  sow=(mjd-44244.d0-week*7)*86400.d0+sod

  RETURN

END SUBROUTINE
