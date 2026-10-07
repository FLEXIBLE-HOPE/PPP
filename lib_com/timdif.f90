!*
REAL(RL) FUNCTION timdif(jd2,sod2,jd1,sod1)
!!
!! purpose   : get time difference
!!
!! parameter :
!!    input  : jd2,sod2 -- time to be differenced
!!             jd1,sod1 -- time difference
!!
!! author    : Geng J.H.
!!
!! created   : 12/10/2006
!*
USE par
IMPLICIT NONE

!*
! The arguments
!!-----------------------
INTEGER(IT) :: jd1,jd2
REAL(RL) :: sod1,sod2

  !*
  ! Start the exectuable code
  !!--------------------------

  timdif=86400.D0*(jd2-jd1)+sod2-sod1

  RETURN

END FUNCTION
