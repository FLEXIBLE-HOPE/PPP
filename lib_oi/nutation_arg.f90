!*
SUBROUTINE nutation_arg(mjd,arg)
!!
!*
USE const
IMPLICIT NONE

!*
! The arguments
!!--------------------
REAL(RL) :: mjd ! TT
REAL(RL) :: arg(1:*)

  !*
  ! The local variables
  !!---------------------------
  INTEGER(IT) :: i,j
  LOGICAL(LG) :: lfirst
  REAL(RL) :: coeff(5,5),t

  ! IERS 2010, P67
  DATA ((coeff(i,j),j=1,5),i=1,5)                                           &
  / 134.96340251d0, 1717915923.2178d0, 31.8792d0, 0.051635d0,-0.00024470d0  &
  , 357.52910918d0,  129596581.0481d0, -0.5532d0, 0.000136d0,-0.00001149d0  &
  ,  93.27209062d0, 1739527262.8478d0,-12.7512d0,-0.001037d0, 0.00000417d0  &
  , 297.85019547d0, 1602961601.2090d0, -6.3706d0, 0.006593d0,-0.00003169d0  &
  , 125.04455501d0,   -6962890.2665d0,  7.4722d0, 0.007702d0,-0.00005939d0/

  DATA lfirst /.TRUE./

  SAVE coeff,lfirst

  !*
  ! Start of executable code
  !!------------------------

  IF (lfirst) THEN
    DO i=1,5
      coeff(i,1)=coeff(i,1)*DEG2RAD
      DO j=2,4
        coeff(i,j)=coeff(i,j)*ARCSEC2RAD
      END DO
    END DO
    lfirst=.FALSE.
  END IF

  ! calculate the angles
  t=(mjd-51544.5d0)/36525.d0
  DO i=1,5
    arg(i)=coeff(i,1)+(coeff(i,2)+(coeff(i,3)+(coeff(i,4)+coeff(i,5)*t)*t)*t)*t
    DO WHILE(arg(i) .LT. -PI)
      arg(i)=arg(i)+PI
    END DO
     DO WHILE(arg(i) .GT. PI)
      arg(i)=arg(i)-PI
    END DO
  END DO


  RETURN

END SUBROUTINE
