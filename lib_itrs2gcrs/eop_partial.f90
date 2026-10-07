!*
SUBROUTINE eop_partial(dt,xsit,dgmat,dxmat,dymat,dpole,dut1)
!!
!*
USE const
IMPLICIT NONE

!*
! The arguments
!!--------------------------
REAL(RL) :: dt,xsit(1:*)
REAL(RL) :: dgmat(3,3),dxmat(3,3)
REAL(RL) :: dymat(3,3),dpole(3,4),dut1(3,2)

  !*
  ! The local variables
  !!------------------------------
  INTEGER(IT) :: i,j

  !*
  ! Start the exectuable code
  !!-----------------------------

  DO i=1, 3
    dpole(i,1)=0.d0
    dpole(i,3)=0.d0
    dut1(i,1)=0.d0
    DO j=1, 3
      dpole(i,1)=dpole(i,1)+dxmat(i,j)*xsit(j)
      dpole(i,3)=dpole(i,3)+dymat(i,j)*xsit(j)
      dut1(i,1)=dut1(i,1)+dgmat(i,j)*xsit(j) !/E_ROTATE*2*PI*1.00273781191135448d0*86400.d0
    END DO
  END DO

  DO i=1, 3
    dpole(i,2)=dpole(i,1)*dt
    dpole(i,4)=dpole(i,3)*dt
    dut1(i,2)=dut1(i,1)*dt
  END DO

  ! From asc to masc
  DO i=1, 3
    DO j=1, 4
      dpole(i,j)=dpole(i,j)*1.d-3
      IF (j .LE. 2) dut1(i,j)=dut1(i,j)*1.d-3
    END DO
  END DO

  RETURN

END SUBROUTINE
