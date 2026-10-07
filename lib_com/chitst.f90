!*
LOGICAL(LG) FUNCTION chitst(flag,ndof,osig,nsig,zone)
!!
!! purpose  : chi-square test
!! parameter:
!!    input : flag -- 0: dual-side test,  H0 nsig^2 == osig^2   H1 nsig^2 != osig^2
!!                    1: right side test, H0 nsig^2 <= osig^2   H1 nsig^2 >  osig^2
!!                   -1: left side test,  H0 nsig^2 >= osig^2   H1 nsig^2 <  osig^2
!!            ndof -- degree of freedom
!!            osig -- reference sigma
!!            nsig -- new sigma
!!            zone -- belief zone 0.95, 0.99 ...
!!    output: chitst -- true, then accept H0, otherwise false, accept H1
!! author   : Geng J
!! created  : Oct. 17, 2007
!! NOTE     : go to page 206 of Chinese book of statistics (Edition of Gaojiao)
!!
!*
USE par
IMPLICIT NONE

!*
! The arguments
!!------------------------
INTEGER(IT) :: flag,ndof
REAL(RL) :: osig,nsig,zone

  !*
  ! The local variables
  !!-------------------------
  INTEGER(IT) :: ierr
  REAL(RL) :: chi2l,chi2r,afa,stat

  !*
  ! Start the exectuable code
  !!------------------------------

  !! chi test
  stat=ndof*(nsig/osig)**2
  chitst=.FALSE.
  IF (flag .EQ. 0) THEN
    afa=(1.d0-zone)/2.d0
    CALL pchi2(ndof,afa,1,chi2r)
    afa=1.d0-afa
    CALL pchi2(ndof,afa,1,chi2l)
    IF (stat.GT.chi2l .AND. stat.LT.chi2r) chitst=.TRUE.
  ELSE IF (flag .GT. 0) THEN
    afa=1.d0-zone
    CALL pchi2(ndof,afa,1,chi2r)
    IF (stat .LT. chi2r) chitst=.TRUE.
  ELSE
    afa=zone
    CALL pchi2(ndof,afa,1,chi2l)
    IF (stat.GT.chi2l) chitst=.TRUE.
  END IF

  RETURN

END FUNCTION
