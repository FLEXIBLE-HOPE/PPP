!*
FUNCTION cprn2sys(prn)
!!
!*
IMPLICIT NONE
CHARACTER(LEN=3), INTENT(IN) :: prn
CHARACTER :: cprn2sys

cprn2sys = prn(1:1)

END FUNCTION cprn2sys
