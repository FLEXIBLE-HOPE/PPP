!*
SUBROUTINE oi_srp_pri(fmodel,mjd,sod,cprn,csvn,blk,PAN,mass,lambda,xsat,xsun,acc)
!!
!*
USE const
USE satellite
IMPLICIT NONE

!*
! The arguments
!!-------------------------
TYPE(SATEPAN) :: PAN
INTEGER(IT) :: mjd
CHARACTER(LEN=*) :: fmodel,cprn,csvn,blk
REAL(RL) :: sod,mass,acc(1:*)
REAL(RL) :: lambda,xsat(1:*),xsun(1:*)

  !*
  ! The local variables
  !!-----------------------------

  INTEGER(IT) :: i,j,k
  REAL(RL) :: factor,cosf,det
  REAL(RL) :: d_unit(3),y_unit(3)
  REAL(RL) :: x_unit(3),z_unit(3)
  REAL(RL) :: sp_unit(3)
  REAL(RL) :: f(3),cost(7)

  !*
  ! The function called
  !!-----------------------------
  REAL(RL) :: dot
  INTEGER(IT) :: pointer_string

  !*
  ! Start the exectuable code
  !!-----------------------------

  ! Vector Sun to satellite
  DO i=1, 3
    d_unit(i)=-(xsat(i)-xsun(i))
  END DO
  CALL unit_vector(3,d_unit,d_unit,det)

  factor=(149597870.6996262d0/det)**2

  CALL rot_scfix2j2000(mjd,sod,cprn,csvn,blk,xsat,xsun,x_unit,y_unit,z_unit)

  factor=-lambda*factor*1367.d0/VEL_LIGHT/mass*1.d-3

  ! take care the code, in the eclipsing season, it is not true
  sp_unit=d_unit

  f=0.d0
  DO i=1, PAN.npan
    ! only for +X, -X, +Y, -Y, +Z, -Z, and solar panels
    SELECT CASE(i)
      CASE(1)
        PAN.normj(1:3,1) = x_unit(1:3)
      CASE(2)
        PAN.normj(1:3,2) =-x_unit(1:3)
      CASE(3)
        PAN.normj(1:3,3) = y_unit(1:3)
      CASE(4)
        PAN.normj(1:3,4) =-y_unit(1:3)
      CASE(5)
        PAN.normj(1:3,5) = z_unit(1:3)
      CASE(6)
        PAN.normj(1:3,6) =-z_unit(1:3)
      CASE(7)
        PAN.normj(1:3,7) = sp_unit(1:3)
    END SELECT

    cost(i)=dot(3,d_unit,PAN.normj(1:3,i))
    IF (DABS(cost(i)) .LT. 1.D-10) cost(i)=0.d0
    IF (cost(i) .GT. 1.D0-1.D-10) cost(i)=1.d0
    IF (cost(i) .LE. 0.d0) CYCLE


    IF (INDEX(fmodel,'PANEL') .NE. 0) THEN
      DO j=1, 3
        ! +X
        IF (i .LE. 6) THEN
          !! Add the Thermal Re-Radiation forces
          f(j)=f(j)+factor*PAN.area(i)*cost(i)* &
               ((PAN.alpha(i)+PAN.delta(i))*(d_unit(j)+2.d0*PAN.normj(j,i)/3.d0)+ &
               2.d0*(PAN.rho(i))*cost(i)*PAN.normj(j,i))
          !! Only Solar Radiation Pressure forces
        ELSE IF (i .EQ. 7) THEN
          f(j)=f(j)+factor*PAN.area(i)*cost(i)* &
                 ((1.d0+PAN.rho(i)+2.d0*PAN.delta(i)/3.d0)*d_unit(j)+2*PAN.rho(i)*(cost(i)*PAN.normj(j,i)-d_unit(j))+ &
                 2.d0*PAN.delta(i)/3.d0*(PAN.normj(j,i)-d_unit(j)))
        END IF
      END DO
    END IF
  END DO

  DO i=1, 3
    acc(i)=acc(i)+f(i)
  END DO

  RETURN

END SUBROUTINE
