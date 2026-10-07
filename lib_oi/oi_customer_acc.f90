!*
SUBROUTINE oi_customer_acc(lpart,cprn,npar,pname,xics,period,dt,acr,cmat,acc0)
!!
!*
USE const
IMPLICIT NONE

!*
! The arguments
!!-----------------------------
LOGICAL(LG) :: lpart
INTEGER(IT) :: npar
CHARACTER(LEN=*) :: cprn,pname(1:*)
REAL(RL) :: xics(1:*),period,dt,acr(3,3)
REAL(RL) :: cmat(1:*),acc0(1:*)

  !*
  ! The local variables
  !!----------------------------
  INTEGER(IT) :: i,j,k,iterm,icomp
  REAL(RL) :: acc(3),term(6),fact

  INTEGER(IT) :: MAXPARLOC
  PARAMETER(MAXPARLOC=18)
  INTEGER(IT) :: ltog(MAXPARLOC)
  CHARACTER(15) :: lpname(MAXPARLOC)
  REAL(RL) :: param(MAXPARLOC)
  DATA lpname &
  /'EMP_Aa    ','EMP_Ba    ','EMP_Sa1   ','EMP_Ca1   ','EMP_Sa2   ','EMP_Ca2   ', &
   'EMP_Ac    ','EMP_Bc    ','EMP_Sc1   ','EMP_Cc1   ','EMP_Sc2   ','EMP_Cc2   ', &
   'EMP_Ar    ','EMP_Br    ','EMP_Sr1   ','EMP_Cr1   ','EMP_Sr2   ','EMP_Cr2   '/

  !*
  ! The function called
  !!----------------------------
  INTEGER(IT) :: pointer_string

  !*
  ! Start the executable code
  !!-----------------------------

  DO i=1, MAXPARLOC
    ltog(i)=pointer_string(npar,pname,lpname(i))
  END DO

  DO i=1, MAXPARLOC
    param(i)=0.d0
    IF (ltog(i) .NE. 0) param(i)=param(i)+xics(ltog(i))
  END DO

  IF (cprn(1:1) .EQ. 'L') THEN
    fact=1.d-9
  ELSE
    fact=1.d-12
  END IF

  ! Constant, linear, one cpr and two cpr
  term(1)=1.d0*fact
  term(2)=dt*fact
  term(3)=DSIN(2*PI/period*dt)*fact
  term(4)=DCOS(2*PI/period*dt)*fact
  term(5)=DSIN(PI/period*dt)*fact
  term(6)=DCOS(PI/period*dt)*fact

  DO i=1, 3
    acc(i)=0.d0
  END DO

  ! CMat and acc
  DO i=1, MAXPARLOC
    IF (ltog(i) .EQ. 0) CYCLE
    icomp=(i-1)/6+1
    iterm=i-(icomp-1)*6
    IF (lpart .EQ. .TRUE.) THEN
      j=(ltog(i)-6-1)*3
      DO k=1, 3
        cmat(j+k)=term(iterm)*acr(k,icomp)
      END DO
    END IF
    IF (param(i) .EQ. 0.d0) CYCLE
    DO k=1,3
      acc(k)=acc(k)+term(iterm)*param(i)*acr(k,icomp)
    END DO
  END DO

  DO i=1, 3
    acc0(i)=acc0(i)+acc(i)
  END DO

  RETURN

END SUBROUTINE
