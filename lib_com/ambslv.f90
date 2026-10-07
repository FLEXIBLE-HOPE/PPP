!*
SUBROUTINE ambslv(ncad,q22,bias,disall)
!!
!! purpose  : resolve ambiguity
!! parameter:
!!    input : ncad   -- # of candidate ambiguities
!!            q22    -- cofactor matrix
!!    output: bias   -- float / fixed ambiguity estimates
!!            disall -- norm of optimum & suboptimum solutions
!! author   : Geng J
!! created  : Mar. 16, 2008
!!
!*
USE par
IMPLICIT NONE

!*
! The arguments
!!--------------------
INTEGER(IT) :: ncad
REAL(RL) :: q22(1:*),bias(1:*),disall(1:*)

  !*
  ! The local variables
  !!---------------------------
  REAL(RL) :: dump

  !*
  ! Start the exectuable code
  !!---------------------------

  IF (ncad.GT.1) THEN
    CALL lambda4(ncad,q22,bias,disall)
  ELSE
    dump=bias(1)
    bias(1)=nint(bias(1))*1.d0
    dump=bias(1)-dump
    disall(1)=dump/q22(1)*dump
    dump=1.d0-dabs(dump)
    disall(2)=dump/q22(1)*dump
  END IF

  RETURN

END SUBROUTINE
