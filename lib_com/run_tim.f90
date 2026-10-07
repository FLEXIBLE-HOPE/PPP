!*
CHARACTER(LEN_RUNTIM) FUNCTION run_tim()
!!
!! purpose  : show running time
!! parameter:
!! author   : Geng J
!! created  : Nov. 13, 2007
!!
!*
USE par
IMPLICIT NONE

  !*
  ! The local variables
  !!------------------------
  INTEGER*2 :: date_time(8)
  CHARACTER(LEN=12) :: real_clock(3)

  !*
  ! Start the exectuable code
  !!--------------------------

  CALL DATE_AND_TIME(real_clock(1),real_clock(2),real_clock(3),date_time)
  WRITE(run_tim,'(I4,1h-,I2.2,1h-,I2.2,A1,3(I2.2,1h:),I4.4,A1)') date_time(1:3),' ',date_time(5:8),' '

  RETURN

END FUNCTION
