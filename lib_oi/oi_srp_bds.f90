!*
SUBROUTINE oi_srp_bds(fmodel,mjd,sod,cprn,csvn,blk,PAN,mass,lambda,xsat,xsun,acc)
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

  REAL(RL) :: dAD(6),dR(6)
  REAL(RL) :: dSP,dSB,dY0


  !*
  ! The function called
  !!-----------------------------
  REAL(RL) :: dot
  INTEGER(IT) :: pointer_string

  !*
  ! Start the exectuable code
  !!-----------------------------

    ! params
  DO i=1,6
    dAD(i)=0.0
    dR(i) =0.0
  END DO
  
  dSP=0.0

  IF (TRIM(cprn) .EQ. 'C19' .OR. TRIM(cprn) .EQ. 'C20' .OR. TRIM(cprn) .EQ. 'C21' .OR. TRIM(cprn) .EQ. 'C22' .OR. TRIM(cprn) .EQ. 'C23' .OR. TRIM(cprn) .EQ. 'C24') THEN 
      dSP=1.19844-1-0.08
      dAD(1)=0.36604-0.35
      dAD(5)=0.10886-0.25
      dAD(6)=0.50315-0.35
      dR(1)=0.60151-0.65
      dR(5)=0.83772-0.75
      dR(6)=0.4643-0.65

  ELSE IF (TRIM(cprn) .EQ. 'C32' .OR. TRIM(cprn) .EQ. 'C33') THEN
      dSP=1.26173-1-0.08
      dAD(1)=0.28049-0.35
      dAD(5)=0.27448-0.25
      dAD(6)=0.42807-0.35
      dR(1)=0.48721-0.65
      dR(5)=0.93970-0.75
      dR(6)=0.45008-0.65

  ELSE IF (TRIM(cprn) .EQ. 'C36' .OR. TRIM(cprn) .EQ. 'C37' .OR. TRIM(cprn) .EQ. 'C41' .OR. TRIM(cprn) .EQ. 'C42') THEN
      dSP=1.24193-1-0.08
      dAD(1)=0.31997-0.35
      dAD(5)=0.10274-0.25
      dAD(6)=0.54920-0.35
      dR(1)=0.52869-0.65
      dR(5)=0.90468-0.75
      dR(6)=0.47575-0.65

  ELSE IF (TRIM(cprn) .EQ. 'C45' .OR. TRIM(cprn) .EQ. 'C46') THEN
      dSP=1.31312-1-0.08
      dAD(1)=0.30465-0.35
      dAD(5)=-0.10019-0.25
      dAD(6)=0.85915-0.35
      dR(1)=0.49278-0.65
      dR(5)=0.94301-0.75
      dR(6)=0.49956-0.65

  ELSE IF (TRIM(cprn) .EQ. 'C25' .OR. TRIM(cprn) .EQ. 'C26') THEN
      dSP=1.26530-1-0.08
      dAD(1)=0.16674-0.2
      dAD(5)=0.42831-0.156
      dAD(6)=0.05916-0.156
      dR(1)=0.77478-0.8
      dR(5)=0.70224-0.844
      dR(6)=0.85930-0.844

  ELSE IF (TRIM(cprn) .EQ. 'C27' .OR. TRIM(cprn) .EQ. 'C28' .OR. TRIM(cprn) .EQ. 'C29' .OR. TRIM(cprn) .EQ. 'C30') THEN
      dSP=1.25299-1-0.08
      dAD(1)=0.15530-0.2
      dAD(5)=0.39270-0.156
      dAD(6)=0.07642-0.156
      dR(1)=0.75653-0.8
      dR(5)=0.72644-0.844
      dR(6)=0.81101-0.844

  ELSE IF (TRIM(cprn) .EQ. 'C34' .OR. TRIM(cprn) .EQ. 'C35') THEN
      dSP=1.37444-1-0.08
      dAD(1)=0.10934-0.2
      dAD(5)=0.49362-0.156
      dAD(6)=0.08586-0.156
      dR(1)=0.68945-0.8
      dR(5)=0.71520-0.844
      dR(6)=0.87707-0.844

  ELSE IF (TRIM(cprn) .EQ. 'C43' .OR. TRIM(cprn) .EQ. 'C44') THEN
      dSP=1.28552-1-0.08
      dAD(1)=0.18562-0.2
      dAD(5)=0.36952-0.156
      dAD(6)=-0.00146-0.156
      dR(1)=0.81968-0.8
      dR(5)=0.71194-0.844
      dR(6)=0.80203-0.844

  ELSE IF (TRIM(cprn) .EQ. 'C38' .OR. TRIM(cprn) .EQ. 'C40') THEN
      dSP=1.08596-1-0.08
      dAD(1)=0.39439-0.35
      dAD(5)=0.68650-0.87
      dAD(6)=0.74121-0.87
      dR(1)=0.75801-0.65
      dR(5)=-0.00994-0.13
      dR(6)=0.02738-0.13

  ELSE IF (TRIM(cprn) .EQ. 'C39') THEN
      dSP=1.11435-1-0.08
      dAD(1)=0.38322-0.35
      dAD(5)=0.58084-0.87
      dAD(6)=0.67228-0.87
      dR(1)=0.65716-0.65
      dR(5)=0.17734-0.13
      dR(6)=-0.00516-0.13

  END IF

  ! Vector Sun to satellite
  DO i=1, 3
    d_unit(i)=-(xsat(i)-xsun(i))
  END DO
  CALL unit_vector(3,d_unit,d_unit,det)

  factor=(149597870.691d0/det)**2

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
               ((PAN.alpha(i)+PAN.delta(i)+dAD(i))*(d_unit(j)+2.d0*PAN.normj(j,i)/3.d0)+ &
               2.d0*(PAN.rho(i)+dR(i))*cost(i)*PAN.normj(j,i))
          !! Only Solar Radiation Pressure forces
          !f(j)=f(j)+factor*PAN.area(i)*cost(i)* &
          !     (d_unit(j)+2.d0/3.d0*(PAN.delta(i))*PAN.normj(j,i)+(2*cost(i)*PAN.normj(j,i)-d_unit(j))*(PAN.rho(i)))
        ELSE IF (i .EQ. 7) THEN
          f(j)=f(j)+factor*PAN.area(i)*cost(i)* &
                 ((1.d0+PAN.rho(i)+2.d0*PAN.delta(i)/3.d0+dSP)*d_unit(j)+2*PAN.rho(i)*(cost(i)*PAN.normj(j,i)-d_unit(j))+ &
                 2.d0*PAN.delta(i)/3.d0*(PAN.normj(j,i)-d_unit(j)))
        END IF
      END DO
    END IF
  END DO

  DO i=1, 3
    acc(i)=acc(i)+f(i)
  END DO
  !IF (sod .EQ. 300.d0) THEN 
  !   write(1200,'(f23.12,1X,A6,1X,A4,7f12.5)') sod,cprn,'islpri',dSP, dAD(1),dAD(5),dAD(6),dR(1),dR(5),dR(6)
  !END IF

  RETURN

END SUBROUTINE

