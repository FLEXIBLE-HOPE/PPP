

REAL(RL) FUNCTION wpress(rh,t)
USE par
IMPLICIT NONE

  REAL(RL) rh,t

  wpress=rh*6.11d0*10.d0**(7.5d0*t/(t+2.373d2))

  RETURN

END FUNCTION
