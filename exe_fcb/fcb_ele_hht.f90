!
!! purpose  : elemental Householder transformation
!! parameter:
!!    input : m -- size of vector a & u
!!            a -- original vector
!!    output: s -- next a(1)
!!            u,beta -- elemental transformation vector
!! author   : Geng J
!! created  : Nov. 11, 2007
!
SUBROUTINE fcb_ele_hht(m,a,s,u,beta)
!!
!*
USE par
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!------------------
INTEGER(IT) :: m
REAL(RL) :: a(1:*),u(1:*),s,beta

  !*
  ! The local variables
  !!-----------------------
  INTEGER(IT) :: i

  !*
  ! Start the exectuable code
  !!--------------------------


  !! compute s
  s=0.d0
  DO i=1,m
    IF (a(i) .EQ. 0.d0) CYCLE
    s=s+a(i)**2
  END DO
  !! use the opposite sign as reference
  s=-dsqrt(s)*dsign(1.d0,a(1))
  IF (s .EQ. 0.d0) THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(ckd_ele_hht): s equal 0.d0'
    CALL exit(1)
  END IF

  !! compute u
  u(1)=a(1)-s
  IF (u(1) .EQ. 0.d0) THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(ckd_ele_hht): u(1) equal 0.d0'
    CALL exit(1)
  END IF

  DO i=2,m
    u(i)=a(i)
  END DO

  !! compute beta
  beta=1.d0/s/u(1)

  RETURN

END SUBROUTINE
