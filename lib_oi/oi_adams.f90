!*
SUBROUTINE oi_adams(str,order,x,f,hh)
!!
!*
USE orbit
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!-----------------------
CHARACTER(LEN=*) :: str
INTEGER(IT) :: order
REAL(RL) :: x(0:*),f(0:*),hh

  !*
  ! The local variables
  !!--------------------------

  INTEGER(IT) :: i
  REAL(RL) :: sumt
  TYPE(ADAMSCOEF) :: adamc

  LOGICAL(LG) :: lfirst
  DATA lfirst /.TRUE./

  SAVE adamc,lfirst

  !*
  ! Start the exectuable code
  !!--------------------

  IF (lfirst .EQ. .TRUE.) THEN
    adamc.n=order
    CALL oi_adams_coef(adamc.n,adamc.pbeta,adamc.cbeta)
    lfirst=.FALSE.
  END IF

  sumt=0.d0
  SELECT CASE(TRIM(str))
    CASE('PRE')
      DO i=0, adamc.n-1
        sumt=sumt+adamc.pbeta(i)*f(adamc.n-1-i)
      END DO
      x(adamc.n)=x(adamc.n-1)+hh*sumt
    CASE('COR')
      DO i=0, adamc.n-1
        sumt=sumt+adamc.cbeta(i)*f(adamc.n-i)
      END DO
      x(adamc.n)=x(adamc.n-1)+hh*sumt
    CASE DEFAULT
      WRITE(ERROR_UNIT,'(A)') '***ERROR(oi_adams): unknown type '//TRIM(str)
      CALL exit(1)
  END SELECT

  RETURN

END SUBROUTINE
