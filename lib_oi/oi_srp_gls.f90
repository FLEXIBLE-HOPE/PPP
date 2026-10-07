!*
SUBROUTINE oi_srp_gls(mjd,sod,cprn,csvn,blk,mass,lambda,xsat,xsun,acc)
!!
!*
USE const
USE satellite
IMPLICIT NONE

!*
! The arguments
!!-------------------------
INTEGER(IT) :: mjd
CHARACTER(LEN=*) :: cprn,csvn,blk
REAL(RL) :: sod,mass,acc(1:*)
REAL(RL) :: lambda,xsat(1:*),xsun(1:*)

  !*
  ! The local variables
  !!-----------------------------

  INTEGER(IT) :: i,j,k
  REAL(RL) :: factor,cosf,det,massx,radiator
  REAL(RL) :: d_unit(3),y_unit(3)
  REAL(RL) :: x_unit(3),z_unit(3)
  REAL(RL) :: sp_unit(3)
  REAL(RL) :: f(3)

  !*
  ! The function called
  !!-----------------------------
  REAL(RL) :: dot
  INTEGER(IT) :: pointer_string

  !*
  ! Start the exectuable code
  !!-----------------------------

  ! the initial values
  SELECT CASE(TRIM(blk))
    CASE('GLONASS-M')
       massx=1415.d0
       radiator=-1.037d-12
    CASE('GLONASS-K1')
       massx=935.d0
       radiator=-0.493d-12
    CASE DEFAULT
       massx=mass
       radiator=0.d0
  END SELECT

  ! Vector Sun to satellite
  DO i=1, 3
    d_unit(i)=-(xsat(i)-xsun(i))
  END DO
  CALL unit_vector(3,d_unit,d_unit,det)

  factor=(149597870.691d0/det)**2
  CALL rot_scfix2j2000(mjd,sod,cprn,csvn,blk,xsat,xsun,x_unit,y_unit,z_unit)
  factor=-lambda*factor*1367.d0/VEL_LIGHT/massx*1.d-3

  ! take care the code, in the eclipsing season, it is not true
  sp_unit=d_unit

  f=0.d0

  ! +X surface
  cosf=dot(3,d_unit,x_unit)
  IF (DABS(cosf) .LT. 1.D-10) cosf=0.d0
  IF (cosf .GT. 1.D0-1.D-10) cosf=1.d0
  IF (cosf .GT. 0.d0) THEN
    DO i=1, 3
      SELECT CASE(TRIM(blk))
        CASE('GLONASS-M')
          f(i)=f(i)+factor*cosf*4.53 * &
             (0.866d0*(d_unit(i)+(PI/6.d0*0.728 +2.d0/3.d0*(1.d0-0.728d0))*x_unit(i)) + &
             (4.d0/3.d0*0.728 + 2.d0*(1.d0-0.728))*0.022*cosf*x_unit(i))
        CASE('GLONASS-K1')
          f(i)=f(i)+factor*cosf*2.21 * &
!             (0.951d0*(d_unit(i)+PI/6.d0*x_unit(i)) - &
!             4.d0/3.d0*0.115*cosf*x_unit(i))
             (0.951d0*(d_unit(i)+2.d0/3.d0*x_unit(i)) - &
             2.d0*0.115*cosf*x_unit(i))
        CASE DEFAULT
          f(j)=0.d0
      END SELECT
    END DO
  END IF

  ! -X
  cosf=dot(3,d_unit,-x_unit)
  IF (DABS(cosf) .LT. 1.D-10) cosf=0.d0
  IF (cosf .GT. 1.D0-1.D-10) cosf=1.d0
  IF (cosf .GT. 0.d0) THEN
    DO i=1, 3
      SELECT CASE(TRIM(blk))
        CASE('GLONASS-M')
          f(i)=f(i)+factor*cosf*4.53 * &
             (0.866d0*(d_unit(i)-(PI/6.d0*0.728 +2.d0/3.d0*(1.d0-0.728d0))*x_unit(i)) - &
             (4.d0/3.d0*0.728 + 2.d0*(1.d0-0.728))*0.022*cosf*x_unit(i))-radiator*x_unit(i)
        CASE('GLONASS-K1')
          f(i)=f(i)+factor*cosf*2.21 * &
!            (0.951d0*(d_unit(i)-PI/6.d0*x_unit(i)) + &
!             4.d0/3.d0*0.115*cosf*x_unit(i))-radiator*x_unit(i)
             (0.951d0*(d_unit(i)-2.d0/3.d0*x_unit(i)) + &
             2.d0*0.115*cosf*x_unit(i))-radiator*x_unit(i)
        CASE DEFAULT
          f(j)=0.d0
      END SELECT
    END DO
  END IF


  ! +Z
  cosf=dot(3,d_unit,z_unit)
  IF (DABS(cosf) .LT. 1.D-10) cosf=0.d0
  IF (cosf .GT. 1.D0-1.D-10) cosf=1.d0
  IF (cosf .GT. 0.d0) THEN
    DO i=1, 3
      SELECT CASE(TRIM(blk))
        CASE('GLONASS-M')
          f(i)=f(i)+factor*cosf*3.40 * &
             (0.479d0*(d_unit(i)+2.d0/3.d0*z_unit(i)) - &
             2.d0*0.169*cosf*z_unit(i))
        CASE('GLONASS-K1')
          f(i)=f(i)+factor*cosf*1.730 * &
!             (0.547d0*(d_unit(i)+PI/6.d0*z_unit(i)) + &
!             4.d0/3.d0*0.217*cosf*z_unit(i))
             (0.547d0*(d_unit(i)+2.d0/3.d0*z_unit(i)) + &
             2.d0*0.217*cosf*z_unit(i))
        CASE DEFAULT
          f(j)=0.d0
      END SELECT
    END DO
  END IF

  ! -Z
  cosf=dot(3,d_unit,-z_unit)
  IF (DABS(cosf) .LT. 1.D-10) cosf=0.d0
  IF (cosf .GT. 1.D0-1.D-10) cosf=1.d0
  IF (cosf .GT. 0.d0) THEN
    DO i=1, 3
      SELECT CASE(TRIM(blk))
        CASE('GLONASS-M')
          f(i)=f(i)+factor*cosf*3.40 * &
             (0.584d0*(d_unit(i)-2.d0/3.d0*z_unit(i)) + &
             2.d0*0.215*cosf*z_unit(i))
        CASE('GLONASS-K1')
          f(i)=f(i)+factor*cosf*1.730 * &
!             (0.533d0*(d_unit(i)-PI/6.d0*z_unit(i)) - &
!             4.d0/3.d0*0.196*cosf*z_unit(i))
             (0.533d0*(d_unit(i)-2.d0/3.d0*z_unit(i)) - &
             2.d0*0.196*cosf*z_unit(i))
        CASE DEFAULT
          f(j)=0.d0
      END SELECT
    END DO
  END IF

  ! solar panel
  cosf=dot(3,d_unit,sp_unit)
  IF (DABS(cosf) .LT. 1.D-10) cosf=0.d0
  IF (cosf .GT. 1.D0-1.D-10) cosf=1.d0
  IF (cosf .GT. 0.d0) THEN
    DO i=1, 3
      SELECT CASE(TRIM(blk))
        CASE('GLONASS-M')
          f(i)=f(i)+factor*30.85*cosf* &
               (0.805*d_unit(i)+2.d0*(0.035/3.d0+0.239*cosf)*sp_unit(i))
        CASE('GLONASS-K1')
          f(i)=f(i)+factor*16.96*cosf* &
               (0.805*d_unit(i)+2.d0*(0.035/3.d0+0.124*cosf)*sp_unit(i))
        CASE DEFAULT
          f(i)=0.d0
      END SELECT
    END DO
  END IF

  DO i=1, 3
    acc(i)=acc(i)+f(i)
  END DO

  RETURN

END SUBROUTINE
