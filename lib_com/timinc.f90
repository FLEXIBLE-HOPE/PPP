!*
SUBROUTINE timinc(jd,sec,delt,jd1,sec1)
!!
!*
USE par
IMPLICIT NONE

!*
! The arguments
!!-----------------------
INTEGER(IT) :: jd,inc,jd1
REAL(RL) :: sec,delt,sec1

  !*
  ! Start the exectuable code
  !!---------------------------

  sec1 = sec + delt
  inc = INT(sec1/86400.D0)
  sec1 = sec1-inc*86400.D0
  jd1 = jd + inc
  IF (sec1 .GE. 0) RETURN
  jd1 = jd1-1
  sec1 = 86400.D0+sec1

  RETURN

END SUBROUTINE
