!*
SUBROUTINE fund_arg_nutation(rmjd_tdt,arg)
!!
!! purpose   : fundamental arguments of nutation
!! 
!! parameters: 
!!             rmjd_tdt -- mjd in TDT system
!!             arg      -- five nutation arguments
!!
!*
USE const
IMPLICIT NONE

!*
! The arguments
!!--------------------
REAL(RL) rmjd_tdt,arg(1:*)

  !*
  ! The local variables
  !!---------------------------
  LOGICAL(LG) first
  INTEGER(IT) i,j
  REAL(RL) coef(5,5),ttc

  !! IERS resolution 2003
  DATA coef &
  / 134.96340251d0, 1717915923.2178d0, 31.8792d0, 0.051635d0,-0.00024470d0  &
  , 357.52910918d0,  129596581.0481d0, -0.5532d0, 0.000136d0,-0.00001149d0  &
  ,  93.27209062d0, 1739527262.8478d0,-12.7512d0,-0.001037d0, 0.00000417d0  &
  , 297.85019547d0, 1602961601.2090d0, -6.3706d0, 0.006593d0,-0.00003169d0  &
  , 125.04455501d0,   -6962890.2665d0,  7.4722d0, 0.007702d0,-0.00005939d0/


  DATA first/.TRUE./
  SAVE coef,first

  !*
  ! Start of executable code
  !!------------------------

  IF (first) THEN
    first=.FALSE.
    DO i=1,5
      coef(1,i)=coef(1,i)*DEG2RAD
      DO j=2,5
        coef(j,i)=coef(j,i)*ARCSEC2RAD
      END DO
    END DO
  END IF

  !! calculate the angles
  ttc=(rmjd_tdt-51544.5d0)/36525.d0
  DO i=1,5
    arg(i)=coef(1,i)+(coef(2,i)+(coef(3,i)+(coef(4,i)+coef(5,i)*ttc)*ttc)*ttc)*ttc
  END DO

  RETURN

END SUBROUTINE
