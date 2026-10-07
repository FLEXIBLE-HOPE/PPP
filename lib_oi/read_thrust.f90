!*
SUBROUTINE read_thrust(t0, t1, sat_type, tman)
!!
!*
IMPLICIT NONE
REAL*8, INTENT(IN) :: t0, t1
CHARACTER(LEN=*), INTENT(IN) :: sat_type
REAL*8, INTENT(OUT) :: tman(2, 6)

tman = 0.d0

END SUBROUTINE read_thrust
