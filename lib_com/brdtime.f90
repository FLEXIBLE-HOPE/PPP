!*
SUBROUTINE brdtime(cprn,mjd,sod)
!!
!*
USE par
IMPLICIT NONE

CHARACTER(LEN_PRN) :: cprn
INTEGER(IT) :: mjd, i
REAL(RL) :: sod, dump(2)

  !*
  ! Start the exectuable code
  !!---------------------------

  SELECT CASE(cprn(1:1))
    CASE('C','B')
      CALL timinc(mjd,sod,14.d0,mjd,sod)
    CASE('R')
      CALL utctai(mjd,sod,mjd,sod)
      CALL timinc(mjd,sod,-19.d0,mjd,sod)
    CASE DEFAULT
  END SELECT

  RETURN

END SUBROUTINE
