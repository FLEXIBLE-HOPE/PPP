!*
SUBROUTINE eci2orb(cprn,xsat,period)
!!
!*
USE const
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!-------------------
REAL(RL) :: xsat(1:*)
REAL(RL) :: period
CHARACTER(LEN_PRN) :: cprn

  !*
  ! The local variables
  !!--------------------------
  REAL(RL) :: r, v, a

  !*
  ! Start the exectuable code
  !!---------------------------

  !! km to m
  r=DSQRT(xsat(1)**2+xsat(2)**2+xsat(3)**2)*1.d3
  v=DSQRT(xsat(4)**2+xsat(5)**2+xsat(6)**2)*1.d3

  a=2.d0/r-v**2/GME

  IF (a .LE. 0.d0) THEN
    WRITE(OUTPUT_UNIT,'(A)') '***WARNNING(eci2orb): the orbit period could not be less than 0 '//TRIM(cprn)
    CALL exit(1)
  END IF

  !! radian/arc
  period=2*PI/(a*DSQRT(a*GME))

  RETURN

END SUBROUTINE
