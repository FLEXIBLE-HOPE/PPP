!*
SUBROUTINE oi_srp_gen(lpart,PAN,mass,lambda,npar,pname,xics,xsat,xsun,acc,cmat)
!!
!*
USE const
USE satellite
IMPLICIT NONE

!*
! The arguments
!!-------------------------
TYPE(SATEPAN) :: PAN
LOGICAL(LG) :: lpart
INTEGER(IT) :: npar
CHARACTER(LEN=*) :: pname(1:*)
REAL(RL) :: mass,xics(1:*),acc(1:*)
REAL(RL) :: lambda,cmat(1:*),xsat(1:*),xsun(1:*)

  !*
  ! The local variables
  !!-----------------------------

  INTEGER(IT) :: i,j,k

  REAL(RL) :: s_unit(3),f(3),cost
  REAL(RL) :: lsun2sat,dist_factor

  INTEGER(IT), PARAMETER :: MAXPARLOC=1
  INTEGER(IT) :: ltog(MAXPARLOC)
  CHARACTER(LEN_ORBPAR) :: lpname(MAXPARLOC)
  REAL(RL) :: param(MAXPARLOC)
  DATA lpname / 'SR_scale  '/

  !*
  ! The function called
  !!-----------------------------
  REAL(RL) :: dot
  INTEGER(IT) :: pointer_string

  !*
  ! Start the exectuable code
  !!-----------------------------

  DO i=1, MAXPARLOC
    ltog(i)=pointer_string(npar,pname,lpname(i))
  END DO

  DO i=1, MAXPARLOC
    param(i)=1.d0
    IF (ltog(i) .NE. 0) param(i)=param(i)+xics(ltog(i))
  END DO

  ! vector satellite to sun
  DO i=1,3
    s_unit(i)=xsun(i)-xsat(i)
  END DO
  CALL unit_vector(3,s_unit,s_unit,lsun2sat)
  dist_factor=(149597870.6996262d0/lsun2sat)**2

  f=0.d0

  DO i=1, PAN.npan
    cost=dot(3,s_unit,PAN.normj(1,i))
    IF (cost .LE. 0.d0) CYCLE
    DO j=1, 3
      f(j)=f(j)-PAN.area(i)*cost*(2.0d0*(PAN.delta(i)/3.0d0+PAN.rho(i)*cost)*PAN.normj(j,i)+(1.d0-PAN.rho(i))*s_unit(j))
    END DO
    ! for satellite but, the thermal radiation should be considered
    IF (INDEX(PAN.name(i),'solar') .NE. 0) CYCLE
    DO j=1, 3
      f(j)=f(j)-PAN.area(i)*cost*2.d0/3.d0*PAN.alpha(i)*PAN%normj(j,i)
    END DO
  END DO

  DO i=1,3
    f(i)=lambda*f(i)/mass*1367.d0*dist_factor/VEL_LIGHT*1.d-3
    acc(i)=acc(i)+param(1)*f(i)
  END DO

  IF (.NOT. lpart) RETURN

  DO i=1, MAXPARLOC
    IF (ltog(i) .EQ. 0) CYCLE
    j=(ltog(i)-6-1)*3
    DO k=1, 3
      SELECT CASE(i)
        CASE(1)
          cmat(j+k)=f(k)
        CASE DEFAULT
      END SELECT
    END DO
  END DO

  RETURN

END SUBROUTINE
