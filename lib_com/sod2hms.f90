!*
SUBROUTINE sod2hms(sod,ih,im,sec)
!!
!! purpose   : seconds of day to hour minute and second
!!
!! created by: Maorong Ge
!*
USE par
IMPLICIT NONE

!*
! The arguments
!!------------------------
INTEGER(IT) :: ih,im
REAL(RL) :: sod,sec

  !*
  ! Start the exectuable code
  !!--------------------------

  ih=INT(sod/3600.d0)
  im=INT((sod-ih*3600.d0)/60.d0)
  sec=sod-ih*3600.d0-im*60.d0

  RETURN

END SUBROUTINE
