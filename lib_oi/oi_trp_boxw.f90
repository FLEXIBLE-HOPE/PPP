!*
SUBROUTINE oi_trp_boxw(lpart,mjd,sod,cprn,csvn,blk,PAN,mass,lambda,npar,pname,xics,xsat,xsun,acc,cmat)
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
INTEGER(IT) :: npar,mjd
CHARACTER(LEN=*) :: cprn,csvn,blk,pname(1:*)
REAL(RL) :: sod,mass,xics(1:*),acc(1:*),amat(3,3)
REAL(RL) :: lambda,cmat(1:*),xsat(1:*),xsun(1:*)

  !*
  ! The local variables
  !!-----------------------------

  INTEGER(IT) :: i,j,k
  REAL(RL) :: factor,cosf,det
  REAL(RL) :: d_unit(3),y_unit(3),b_unit(3)
  REAL(RL) :: x_unit(3),s_unit(3),z_unit(3)
  REAL(RL) :: sp_unit(3)
  REAL(RL) :: t(8),f(3),cost(7)

  REAL(RL) :: f_dyb(3),f_xyz(3),f_acr(3)
  REAL(RL) :: a_unit(3),c_unit(3),r_unit(3)

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

  factor=(149597870.691d0/det)**2*1367.d0

  CALL rot_scfix2j2000(mjd,sod,cprn,csvn,blk,xsat,xsun,x_unit,y_unit,z_unit)

  !CALL cross(d_unit,y_unit,b_unit)
  !CALL unit_vector(3,b_unit,b_unit,det)

  !CALL cross(y_unit,b_unit,sp_unit)
  !CALL unit_vector(3,sp_unit,sp_unit,det)

  !DO i=1 ,3
  !  z_unit(i)=-xsat(i)
  !END DO
  !CALL unit_vector(3,z_unit,z_unit,det)
  !CALL cross(z_unit,d_unit,y_unit)
  !CALL unit_vector(3,y_unit,y_unit,det)
  !CALL cross(y_unit,z_unit,x_unit)
  !CALL unit_vector(3,x_unit,x_unit,det)
  !sp_unit=d_unit

  sp_unit=d_unit

  f=0.d0
  DO i=1, PAN.npan
    ! only for +X, -X, +Y, -Y, +Z, -Z, and solar panels
    SELECT CASE(i)
      CASE(1)
        PAN.norm(1:3,1) = x_unit(1:3)
      CASE(2)
        PAN.norm(1:3,2) =-x_unit(1:3)
      CASE(3)
        PAN.norm(1:3,3) = y_unit(1:3)
      CASE(4)
        PAN.norm(1:3,4) =-y_unit(1:3)
      CASE(5)
        PAN.norm(1:3,5) = z_unit(1:3)
      CASE(6)
        PAN.norm(1:3,6) =-z_unit(1:3)
      CASE(7)
        PAN.norm(1:3,7) = sp_unit(1:3)
    END SELECT

    cost(i)=dot(3,d_unit,PAN.norm(1:3,i))
    IF (DABS(cost(i)) .LT. 1.D-10) cost(i)=0.d0
    IF (cost(i) .GT. 1.D0-1.D-10) cost(i)=1.d0
    IF (cost(i) .LT. 0.d0) cost(i)=0.d0

    IF (i .EQ. 7) THEN
      CALL solar(factor*cost(i),PAN.area(i)/2.d0,PAN.alpha(i),t(7:8))
    ELSE
      CALL bus(factor*cost(i),PAN.alpha(i),t(i))
    END IF

    DO j=1, 3
      IF (i .NE. 7) THEN
        f(j)=f(j)-2.d0/3.d0*PAN.area(i)*5.670373d-8*t(i)**4/VEL_LIGHT/mass*0.72
      ELSE
        f(j)=f(j)-2.d0/3.d0*PAN.area(i)*5.670373d-8*(t(7)**4-t(8)**4)/VEL_LIGHT/mass*0.86
      !   f(j)=f(j)-factor*PAN.area(i)/VEL_LIGHT/mass*cost(i)*2.d0/3.d0*PAN.alpha(i)*PAN.norm(j,i)
      END IF
    END DO
  END DO

  DO i=1, 3
    acc(i)=acc(i)+f(i)*1.d-3
  END DO

    !if (cprn(1:3) .eq. 'C06') then
    !DO i=1, 3
    !  r_unit(i)=xsat(i)
    !  a_unit(i)=xsat(3+i)
    !END DO

    !CALL unit_vector(3,r_unit,r_unit,det)
    !CALL unit_vector(3,a_unit,a_unit,det)
    !CALL cross(r_unit,a_unit,c_unit)
    !CALL unit_vector(3,c_unit,c_unit,det)
    !CALL cross(c_unit,r_unit,a_unit)
    !CALL unit_vector(3,a_unit,a_unit,det)
    !f_acr(1)=dot(3,f,a_unit)
    !f_acr(2)=dot(3,f,c_unit)
    !f_acr(3)=dot(3,f,r_unit)

    !f_xyz(1)=dot(3,f,x_unit)
    !f_xyz(2)=dot(3,f,y_unit)
    !f_xyz(3)=dot(3,f,z_unit)

    !CALL cross(d_unit,z_unit,y_unit)
    !CALL unit_vector(3,y_unit,y_unit,det)
    !CALL cross(d_unit,y_unit,b_unit)
    !CALL unit_vector(3,b_unit,b_unit,det)
    !f_dyb(1)=dot(3,f,d_unit)
    !f_dyb(2)=dot(3,f,y_unit)
    !f_dyb(3)=dot(3,f,b_unit)

    !write(2000,'(f23.12,1X,A3,1X,20E23.12)') mjd+sod/86400.d0,cprn,f,f_acr,f_dyb,f_xyz,DSQRT(f(1)**2+f(2)**2+f(3)**2),cost

    !write(3000,'(f23.12,1X,A3,1X,20E23.12)') mjd+sod/86400.d0,cprn,t
  !end if


  IF (lpart .EQ. .FALSE.) RETURN

  RETURN

END SUBROUTINE


SUBROUTINE solar(w,are,alp,t)
!!
USE const
IMPLICIT NONE

!*
! The arguments
!!------------------
REAL(RL) ::w,are,alp,t(2)

  !*
  ! The local variables
  !!-----------------------------
  REAL(RL) :: eps,dlt,rto
  REAL(RL) :: qec,e,x,dx,f,df

  DATA eps,dlt,rto,qec /0.86d0,5.670373d-8,0.016055929d0,90.0d0/

  !*
  ! Start the exectuable code
  !!-----------------------------

  x=272.d0
  dx=1.d0

  DO WHILE (DABS(dx) .GT. 1.d-5)
    f=eps*dlt*((x+rto*eps*dlt*are*x**4)**4)+eps*dlt*x**4+qec-alp*w
    df=4.d0*eps*dlt*((x+rto*eps*dlt*are*x**4)**3)*(1.d0+4.d0*rto*eps*dlt*are*x**3)+4.d0*eps*dlt*x**3
    dx=f/df
    x=x-dx
  END DO

  t(1)=x
  t(2)=x+rto*eps*dlt*are*x**4

  RETURN

END SUBROUTINE


SUBROUTINE bus(w,alp,t)
!!
USE const
IMPLICIT NONE

!*
! The arguments
!!----------------------
REAL(RL) :: w,alp,t

  !*
  ! The local variables
  !!-------------------------
  REAL(RL) :: eml,eff,dlt

  DATA eml,eff,dlt /0.72d0,0.02d0,5.670373d-8/

  !*
  ! Start the exectuable code
  !!--------------------------

  t=DSQRT(DSQRT((alp*w+eff*dlt*298.15d0**4)/(dlt*(eml+eff))))

  RETURN

END SUBROUTINE
