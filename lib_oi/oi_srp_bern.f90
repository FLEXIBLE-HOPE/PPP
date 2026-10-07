!*
!! to be test for beidou, and established beidou srp
SUBROUTINE oi_srp_bern(fmodel,lpart,mjd,sod,cprn,csvn,blk,PAN,mass,lambda,npar,pname,xics,xsat,xsun,acc,cmat)
!!
!*
USE par
USE satellite
IMPLICIT NONE

!*
! The arguments
!!-------------------------
TYPE(SATEPAN) :: PAN
LOGICAL(LG) :: lpart
INTEGER(IT) :: npar,mjd
CHARACTER(LEN=*) :: fmodel,cprn,csvn,blk,pname(1:*)
REAL(RL) :: sod,mass,xics(1:*),acc(1:*)
REAL(RL) :: lambda,cmat(1:*),xsat(1:*),xsun(1:*)

  !*
  ! The local variables
  !!------------------------
  INTEGER(IT) :: i,j,k

  REAL(RL) :: factor,det,eps
  REAL(RL) :: xsat_unit(3),node_unit(3)
  REAL(RL) :: d_unit(3),y_unit(3),b_unit(3)
  REAL(RL) :: x_unit(3),s_unit(3),z_unit(3)
  REAL(RL) :: nop_unit(3),xsun_unit(3)
  REAL(RL) :: z_inertial(3),d0_const
  REAL(RL) :: fx,fz,u,u0,f(3),fp(3),cost
  REAL(RL) :: d_d,d_y,d_b,d_x,d_xp1,d_xp3,d_zp
  
  REAL(RL) :: beta,mu

  INTEGER(IT), PARAMETER :: MAXPARLOC=16
  INTEGER(IT) :: ltog(MAXPARLOC)
  CHARACTER(LEN_ORBPAR) :: lpname(MAXPARLOC)
  REAL(RL) :: param(MAXPARLOC)
  DATA lpname &
  /'Kd_BERN   ','Ky_BERN   ','Kb_BERN   ', &
   'Kxp1_BERN ','Kxp3_BERN ','Kzp_BERN  ', &
   'Kdc1_BERN ','Kds1_BERN ','Kyc1_BERN ', &
   'Kys1_BERN ','Kbc1_BERN ','Kbs1_BERN ', &
   'Kdc2_BERN ','Kds2_BERN ','Kdc4_BERN ', &
   'Kds4_BERN '/

  DATA z_inertial /0.d0,0.d0,1.d0/

  REAL(RL) :: a_unit(3),c_unit(3),r_unit(3),f_acr(3),f_xyz(3),f_dyb(3)

  !real(rl) sods
  !integer(it) mjds
  !character(len_prn) cprns
  !logical(lg) lfirst,lacc
  !save sods,mjds,lfirst

  !*
  ! The function used
  !!----------------------------
  REAL(RL) :: dot
  INTEGER(IT) :: pointer_string

  !*
  ! Start the exectuable code
  !!----------------------------

  DO i=1, MAXPARLOC
    ltog(i)=pointer_string(npar,pname,lpname(i))
  END DO

  DO i=1, MAXPARLOC
    param(i)=0.d0
    IF (ltog(i) .NE. 0) param(i)=param(i)+xics(ltog(i))
  END DO

  ! Z axis of satellite
  CALL unit_vector(3,xsat,xsat_unit,det)
  DO i=1, 3
    z_unit(i)=-xsat_unit(i)
  END DO

  ! Vector Sun to satellite
  DO i=1, 3
    s_unit(i)=xsat(i)-xsun(i)
  END DO
  CALL unit_vector(3,s_unit,s_unit,det)
  DO i=1, 3
    d_unit(i)=-s_unit(i)
  END DO
  eps=DACOS(dot(3, d_unit, z_unit))

  factor=(149597870.691d0/det)**2
  
  CALL betau(xsat,xsun,beta,u)
  !! u0 should be defined
  u0=0.d0
 
  !! To get the Y-axis of satellite reference
  !! BeiDou IGSO/MEO/GEO in ON should be tested
  !CALL rot_scfix2j2000(mjd,sod,cprn,blk,xsat,xsun,x_unit,y_unit,z_unit)
  CALL rot_scfix2j2000(mjd,sod,cprn,csvn,blk,xsat,xsun,x_unit,y_unit,z_unit)

  !IF (TRIM(cprn) .EQ. 'C01') THEN
  IF (TRIM(blk) .EQ. 'BEIDOU-2G') THEN
    DO i=1, 3
      z_unit(i)=-xsat_unit(i)
    END DO

    CALL cross(z_unit,d_unit,y_unit)
    CALL unit_vector(3,y_unit,y_unit,det)

    CALL cross(y_unit,z_unit,x_unit)
    CALL unit_vector(3,x_unit,x_unit,det)
  END IF
!????????????????????
!  IF (cprn(1:1).EQ.'G' .OR. cprn(1:1).EQ.'R') THEN
!    DO i=1, 3
!      z_unit(i)=-xsat_unit(i)
!    END DO
!
!    CALL cross(z_unit,d_unit,y_unit)
!    CALL unit_vector(3,y_unit,y_unit,det)
!
!    CALL cross(y_unit,z_unit,x_unit)
!    CALL unit_vector(3,x_unit,x_unit,det)
!  END IF

  !! With Y axis to get D and B
  CALL cross(d_unit,y_unit,b_unit)
  CALL unit_vector(3,b_unit,b_unit,det)
  CALL cross(y_unit,b_unit,d_unit)
  CALL unit_vector(3,d_unit,d_unit,det)

  !! Unit scale factor
  SELECT CASE(TRIM(blk))
    CASE('BLOCK I')
      d0_const=4.54d-5
    CASE('BLOCK II')
      d0_const=8.695d-5
    CASE('BLOCK IIA')
      d0_const=8.958d-5
    CASE('BLOCK IIR-A')
      d0_const=11.17d-5
    CASE('BLOCK IIR-B')
      d0_const=11.17d-5
    CASE('BLOCK IIR-M')
      d0_const=10.88d-5
    CASE('BLOCK IIF')
      d0_const=16.70d-5
    CASE('BLOCK IIIA')
      d0_const=17.50d-5
    CASE('GLONASS')
      d0_const=mass
    CASE('GLONASS-M')
      d0_const=21.03d-5
    CASE('GLONASS-K1')
      d0_const=7.95d-5
    CASE('GALILEO-0A')
      d0_const=mass
    CASE('GALILEO-0B')
      d0_const=mass
    CASE('GALILEO-1')
      d0_const=12.5d-5
    CASE('GALILEO-2')
      d0_const=12.5d-5
    CASE('BEIDOU-2G')
      d0_const=17.76d-5
    CASE('BEIDOU-2I')
      d0_const=17.18d-5
    CASE('BEIDOU-2M')
      d0_const=14.27d-5
    CASE('BEIDOU-3IS-SECM','BEIDOU-3SI-SECM')
      d0_const=mass
    CASE('BEIDOU-3IS-CAST','BEIDOU-3SI-CAST')
      d0_const=27.18d-5
    CASE('BEIDOU-3MS-SECM','BEIDOU-3SM-SECM')
      d0_const=mass
    CASE('BEIDOU-3MS-CAST','BEIDOU-3SM-CAST')
      d0_const=14.32d-5
    CASE('BEIDOU-3G-SECM')
      d0_const=mass
    CASE('BEIDOU-3G-CAST')
      d0_const=mass
    CASE('BEIDOU-3I-SECM')
      d0_const=mass
    CASE('BEIDOU-3I-CAST')
      d0_const=36.98d-5
    CASE('BEIDOU-3M-SECM')
      d0_const=8.59d-5
    CASE('BEIDOU-3M-CAST')
      d0_const=13.79d-5
    CASE('QZSS')
      d0_const=38.04d-5
    CASE('QZSS-2I')
      d0_const=27.04d-5
    CASE('QZSS-2G')
      d0_const=27.04d-5
    CASE('QZSS-2A')
      d0_const=27.04d-5
    CASE('IRNSS-1IGSO')
      d0_const=11.89d-5
    CASE('IRNSS-1GEO')
      d0_const=11.51d-5
    CASE DEFAULT
      d0_const=1367.d0/VEL_LIGHT
  END SELECT
  d0_const=d0_const/mass


  !! The prior model from adjustable coefficients
  d_d=0.d0
  d_y=0.d0
  d_b=0.d0
  d_x=0.d0
  d_xp1=0.d0
  d_xp3=0.d0
  d_zp=0.d0

  cost=dot(3,d_unit,x_unit)
  IF (cost .GT. 0.d0) then
    d_x=cost
  else
    d_x=0.d0
  end if

  IF (ltog( 1) .NE. 0) d_d   = d0_const*param(1)
  IF (ltog( 2) .NE. 0) d_y   = d0_const*param(2)
  IF (ltog( 3) .NE. 0) d_b   = d0_const*param(3)
  IF (ltog( 4) .NE. 0) d_xp1 = param(4)*1.d-6*d_x   ! d_xp1 = d0_const*param(4)*DSIN(u-u0)
  IF (ltog( 5) .NE. 0) d_xp3 = d0_const*param(5)*DSIN(3*u-u0)
  IF (ltog( 6) .NE. 0) d_zp  = d0_const*param(6)*DSIN(u-u0)
  IF (ltog( 7) .NE. 0) d_d = d_d+d0_const*param( 7)*DCOS(u)
  IF (ltog( 8) .NE. 0) d_d = d_d+d0_const*param( 8)*DSIN(u)
  IF (ltog( 9) .NE. 0) d_y = d_y+d0_const*param( 9)*DCOS(u)
  IF (ltog(10) .NE. 0) d_y = d_y+d0_const*param(10)*DSIN(u)
  IF (ltog(11) .NE. 0) d_b = d_b+d0_const*param(11)*DCOS(u)
  IF (ltog(12) .NE. 0) d_b = d_b+d0_const*param(12)*DSIN(u)
  IF (ltog(13) .NE. 0) d_d = d_d+d0_const*param(13)*DCOS(2*u)
  IF (ltog(14) .NE. 0) d_d = d_d+d0_const*param(14)*DSIN(2*u)
  IF (ltog(15) .NE. 0) d_d = d_d+d0_const*param(15)*DCOS(4*u)
  IF (ltog(16) .NE. 0) d_d = d_d+d0_const*param(16)*DSIN(4*u)


  !! The prior solar radiation pressure model
  ! 1) Galileo IOV and FOC satellites
  ! O. Montenbruck, P. Steigenberger, U. Hugentobler (2014)
  ! Enhanced solar radiation pressure modeling for Galileo satellites,
  ! Journal of Geodesy
  ! P. Steigenberger, O. Montenbruck (2016)
  ! Galileo status: orbits, clocks, and positioning
  ! GPS Solutions
  fp=0.d0
  !IF (TRIM(blk).EQ.'BLOCK I' .OR. TRIM(blk).EQ.'BLOCK II' .OR. TRIM(blk).EQ.'BLOCK IIA' .OR. &
  !    TRIM(blk).EQ.'BLOCK IIIA') THEN
  IF (cprn(1:1).EQ.'G' .OR. cprn(1:1).EQ.'R') THEN ! .OR. cprn(1:1).EQ.'J') THEN

    CALL oi_srp_pri(fmodel,mjd,sod,cprn,csvn,blk,PAN,mass,lambda,xsat,xsun,fp)

  !ELSE IF (TRIM(blk).EQ.'GLONASS-K1' .OR. TRIM(blk).EQ.'GLONASS-M') THEN

  !  CALL oi_srp_gls(mjd,sod,cprn,csvn,blk,mass,lambda,xsat,xsun,fp)

  ELSE IF (TRIM(blk).EQ.'GALILEO-1' .OR. TRIM(blk).EQ.'GALILEO-2') THEN
    d_d=d_d-(14.5d-9*(DABS(DCOS(eps))+DSIN(eps)+2.d0/3.d0)+ &
              5.0d-9*(DABS(DCOS(eps))-DSIN(eps)-4.d0/3.d0*DSIN(eps)**2+2.d0/3.d0))+ &
    !! The prior d0 value is taken from Steigenber and Montenbruck (2016)
    ! Steigenber P, Montenbruck O (2016) Galileo status: orbits, clocks, and positioning. GPS Solutions
    ! as the difference between FOC and IOV is less, hence IOV value is used
             87.0d-9
    d_b=d_b-4.d0/3.d0*5.0d-9*(DCOS(eps)*DSIN(eps))

    ! UCL model from Zhen Li
    !CALL uclsrp(-x_unit,-y_unit,z_unit,d_unit,xsat,xsun,lambda,fp)
    !do i=1, 3
    !  fp(i)=fp(i)*1.d-3
    !end do


  !ELSE IF (TRIM(cprn) .EQ. 'C01') THEN
  ! Although the model is developed for C01, but other GEO also show better ODD.
  ELSE IF (TRIM(blk) .EQ. 'BEIDOU-2G') THEN
    !! Wang C. SRP ON mode
    !d_d=d_d+(-0.16d-9*DABS(beta*180.d0/PI)+3.16d-9)*DCOS(u)- &
    !         10.68d-9*DCOS(2.d0*u)- &
    !          1.41d-9*DCOS(4.d0*u)
    !d_y=d_y+1.42d-9*(beta*180.d0/PI)
    !d_b=d_b+(-0.20d-9*DABS(beta*180.d0/PI)+4.42d-9)*DCOS(u)- &
    !          5.93d-9*DCOS(2.d0*u)- &
    !          3.41d-9*DCOS(4.d0*u)

    !IF (DABS(beta) .GE. 8.7d0*PI/180.d0) THEN
    !   d_d=d_d+0.857d-9*(DABS(beta*180.d0/PI)-8.7d0)-113.1d-9
    !ELSE
    !   d_d=d_d-113.1d-9
    !END IF


    !! Wang C. SRP YS mode without shadow
    !d_d=d_d+1.099d-9*DCOS(eps)- &
    !       11.690d-9*DCOS(2.d0*eps)- &
    !        1.460d-9*DCOS(4.d0*eps)
    !d_y=d_y-0.386d-9*(beta*180.d0/PI)*DSIN(u)
    !d_b=d_b+1.27d-9+(-0.448d-9*DABS(beta*180.d0/PI)-4.95d-9)*DCOS(eps)- &
    !          1.25d-9*DCOS(3.d0*eps)

    !IF (DABS(beta) .GE. 8.7d0*PI/180.d0) THEN
    !   d_d=d_d+0.465d-9*(DABS(beta*180.d0/PI)-8.7d0)-113.0d-9
    !ELSE
    !   d_d=d_d-113.0d-9
    !END IF

    !! YS mode with shadow for PRN C01 (SVN C003)
!    d_d=d_d+0.0136d-9*(beta*180.d0/PI)**2-112.1d-9- &
!            0.2200d-9*DCOS(eps)- &
!           10.7000d-9*DCOS(2.d0*eps)- &
!            2.3800d-9*DCOS(4.d0*eps)
!    d_y=d_y-0.4160d-9*(beta*180.d0/PI)*DSIN(u)
!    d_b=d_b+(-0.02167d-9*(beta*180.d0/PI)**2-6.536d-9)*DCOS(eps)- &
!              1.05d-9*DCOS(3.d0*eps)

    !! YS mode with shadow for PRN C01 (SVN GEO08)
    d_d=d_d+0.0146d-9*(beta*180.d0/PI)**2-103.1d-9- &
            1.2200d-9*DCOS(eps)- &
            9.8500d-9*DCOS(2.d0*eps)- &
            2.3600d-9*DCOS(4.d0*eps)
    d_y=d_y-0.520d-9*(beta*180.d0/PI)*DSIN(u)-1.9d-9*DSIN(3.d0*u)
    d_b=d_b+(-0.022d-9*(beta*180.d0/PI)**2-0.8d-9)*DCOS(eps)- &
              2.5d-9*DCOS(3.d0*eps)

!  ELSE IF (INDEX(blk,'BEIDOU-3') .NE. 0) THEN
!    CALL oi_srp_bds(fmodel,mjd,sod,cprn,csvn,blk,PAN,mass,lambda,xsat,xsun,fp)
!!    CALL oi_srp_pri(fmodel,mjd,sod,cprn,csvn,blk,PAN,mass,lambda,xsat,xsun,fp)

   
  ELSE IF (TRIM(blk) .EQ. 'BEIDOU-3M-CAST') THEN

    d_d=d_d-0.00085d-9*(beta*180.d0/PI)**2-140.0d-9+ &
            ! without considering ERP and AT
            !0.5000d-9*DCOS(u)+ &
            ! need to consider ERP and AT
            !0.7000d-9*DCOS(u)+ &
            1.0000d-9*DCOS(u)+ &
            (-0.00039d-9*(beta*180.d0/PI)**2+2.9d-9)*DCOS(2.d0*u)- &
            0.5d-9*DCOS(4.d0*u)
    d_b=d_b+2.3d-9*DCOS(u)+ &
            (0.00058d-9*(beta*180.d0/PI)**2-3.5d-9)*DCOS(3.d0*u)

    ! DLR model
!    d_d=d_d-(11.7d-9*(DABS(DCOS(eps))+DSIN(eps)+2.d0/3.d0)- &
!              5.0d-9*(DABS(DCOS(eps))-DSIN(eps)-4.d0/3.d0*DSIN(eps)**2+2.d0/3.d0)- &
!              0.1d-9*(DCOS(eps)+2.d0/3.d0*DABS(DCOS(eps))*DCOS(eps)))
!    d_b=d_b-(-4.d0/3.d0*5.0d-9*(DCOS(eps)*DSIN(eps))-2.d0/3.d0*0.1d-9*DABS(DCOS(eps))*DSIN(eps))

  ELSE IF (TRIM(blk) .EQ. 'BEIDOU-3M-SECM') THEN

    d_d=d_d+0.0018d-9*(beta*180.d0/PI)**2-75.2d-9- &
            ! without considering ERP and AT
            !2.4000d-9*DCOS(u)+ &  3.0
            ! need to consider ERP and AT
            0.8000d-9*DCOS(u)+ &
            (0.0032d-9*(beta*180.d0/PI)**2-5.2d-9)*DCOS(2.d0*u)- &
            !(0.0032d-9*(beta*180.d0/PI)**2-6.1d-9)*DCOS(2.d0*u)- &
            0.90d-9*DCOS(4.d0*u)
    d_b=d_b-3.5d-9*DCOS(u)- &
    !d_b=d_b-4.5d-9*DCOS(u)- &
            ! without considering ERP and AT
            !(0.0007d-9*(beta*180.d0/PI)**2+1.d-9)*DCOS(3.d0*u)
            ! need to consider ERP and AT
            1.4d-9*DCOS(3.d0*u)

    ! DLR model
!    d_d=d_d-(8.6d-9*(DABS(DCOS(eps))+DSIN(eps)+2.d0/3.d0)+ &
!             1.5d-9*(DABS(DCOS(eps))-DSIN(eps)-4.d0/3.d0*DSIN(eps)**2+2.d0/3.d0)+ &
!             1.4d-9*(DCOS(eps)+2.d0/3.d0*DABS(DCOS(eps))*DCOS(eps)))
!    d_b=d_b-(4.d0/3.d0*1.5d-9*(DCOS(eps)*DSIN(eps))+2.d0/3.d0*1.4d-9*DABS(DCOS(eps))*DSIN(eps))

  ELSE IF (TRIM(cprn) .EQ. 'C06') THEN
    d_xp1=d_xp1-1.5d-9*d_x
  ELSE IF (TRIM(cprn) .EQ. 'C07') THEN
    d_xp1=d_xp1-1.6d-9*d_x
  ELSE IF (TRIM(cprn) .EQ. 'C08') THEN
    d_xp1=d_xp1-1.2d-9*d_x
  ELSE IF (TRIM(cprn) .EQ. 'C09') THEN
    d_xp1=d_xp1-2.3d-9*d_x
  ELSE IF (TRIM(cprn) .EQ. 'C10') THEN
    d_xp1=d_xp1-1.8d-9*d_x
  ELSE IF (TRIM(cprn) .EQ. 'C11') THEN
    d_xp1=d_xp1-1.6d-9*d_x
  ELSE IF (TRIM(cprn) .EQ. 'C12') THEN
    d_xp1=d_xp1-1.5d-9*d_x
  ELSE IF (TRIM(cprn) .EQ. 'C13') THEN
    d_xp1=d_xp1-2.6d-9*d_x
  ELSE IF (TRIM(cprn) .EQ. 'C14') THEN
    d_xp1=d_xp1-0.9d-9*d_x
  
!  ELSE IF (TRIM(cprn) .EQ. 'J01') THEN
  ELSE IF (cprn(1:1) .EQ. 'J') THEN
    ! Montenbruck O, Steigenberger P, Darugna F (2017) Semi-analytical solar radiation pressure modeling for QZS-1
    ! orbit-normal and yaw-steering attitude. Advances in Space Research.
    !IF (DABS(beta*180.0/PI) .GT. 20.0d0)THEN
    !  d_d=d_d-(20.0d-9*(DABS(DCOS(eps))+DSIN(eps)+2.d0/3.d0)- &
    !            7.0d-9*(DABS(DCOS(eps))-DSIN(eps)-4.d0/3.d0*DSIN(eps)**2+2.d0/3.d0)+ &
    !            7.0d-9*DABS(0.5d0*DSIN(2.d0*beta))+ &
    !          112.5d-9)
    !  d_b=d_b+4.d0/3.d0*7.0d-9*(DCOS(eps)*DSIN(eps))
    !ELSE
    !  d_d=d_d-(20.0d-9*((DABS(DCOS(u))+DABS(DSIN(u)))*DCOS(beta)**2+2.d0/3.d0*DCOS(beta))- &
    !            7.0d-9*((DABS(DCOS(u))-DABS(DSIN(u)))*DCOS(beta)**2+2.d0/3.d0*DCOS(beta)*DCOS(2.d0*u))+ &
    !          112.5d-9*DCOS(beta)**2)
    !  d_b=d_b+4.d0/3.d0*7.0d-9*(DCOS(u)*DSIN(u)*DCOS(beta))
    !  d_y=d_y+(20.0d-9*(DABS(DCOS(u))+DABS(DSIN(u)))*0.5d0*DSIN(2.d0*beta)- &
    !            7.0d-9*(DABS(DCOS(u))-DABS(DSIN(u)))*0.5d0*DSIN(2.d0*beta)+ &
    !            7.0d-9*(DABS(DSIN(beta))+2.d0/3.d0)*DSIN(beta)+ &
    !            2.d0*15.d-9*DABS(DSIN(beta))*DSIN(beta)+ & 
    !           70.5d-9*0.5d0*DSIN(2.d0*beta))
    !END IF

    ! Zhao Q, Guo C, Guo J, Liu J, Liu X (2017) An a priori solar radiation pressure model for the QZSS
    ! Michibiki satellite. Journal of Geodesy
    d_d=d_d-26.9d-9*(DABS(DCOS(eps))+DSIN(eps)+2.0d0/3.0d0)+ &
            13.0d-9*(DABS(DCOS(eps))-DSIN(eps)-4.0d0/3.0d0*DSIN(eps)**2+2.0d0/3.0d0)- &
             1.6d-9*(DCOS(eps)+2.0d0/3.0d0*DABS(DCOS(eps))*DCOS(eps))- &
             3.3d-9*(DABS(DCOS(eps))*DCOS(eps)**2+DSIN(eps)**3)- &
        2.d0*3.8d-9*(DABS(DCOS(eps))*DCOS(eps)**2-DSIN(eps)**3)+ &
        2.d0*1.5d-9*DCOS(eps)**3- &
            93.5d-9
    d_b=d_b+4.d0/3.d0*13.0d-9*DCOS(eps)*DSIN(eps)- &
            2.d0/3.d0*1.6d-9*DABS(DCOS(eps))*DSIN(eps)- &
            2.d0*3.3d-9*(DABS(DCOS(eps))-DSIN(eps))*DCOS(eps)*DSIN(eps)- &
            2.d0*3.8d-9*(DABS(DCOS(eps))+DSIN(eps))*DCOS(eps)*DSIN(eps)+ &
            2.d0*1.5d-9*DCOS(eps)**2*DSIN(eps)

    IF (cprn(1:3) .EQ. 'J01') THEN
      IF (DABS(beta*180.0/PI).LE.20.0d0)THEN
        d_d=d_d+128.5d-9*beta**2
        d_b=d_b-3.0d-9
        d_y=d_y+109.0d-9*beta
      END IF
    END IF

!  ELSE IF (TRIM(blk).EQ.'QZSS-2I' .OR. TRIM(blk).EQ.'QZSS-2G') THEN
!
!    CALL oi_srp_pri(fmodel,mjd,sod,cprn,csvn,blk,PAN,mass,lambda,xsat,xsun,fp)

  END IF

  f=0.d0
  DO i=1, 3
    f(i)=lambda*d_d*d_unit(i)+d_y*y_unit(i)+d_b*b_unit(i)+ &
           (d_xp1+d_xp3)*x_unit(i)+d_zp*z_unit(i)
  END DO

  DO i=1, 3
    f(i)=f(i)*factor
    acc(i)=acc(i)+f(i)*1.d-3+fp(i)
  END DO


  IF (.NOT. lpart) RETURN

  !! Acceleration in m/s**2
  DO i=1, MAXPARLOC
    IF (ltog(i) .EQ. 0) CYCLE
    j=(ltog(i)-6 -1)*3
    DO k=1, 3
      IF(i .EQ. 1) THEN
        cmat(j+k)=factor*1.d-3*d0_const*lambda*d_unit(k)
      ELSE IF (i .EQ. 2) THEN
        cmat(j+k)=factor*1.d-3*d0_const*y_unit(k)
      ELSE IF (i .EQ. 3) THEN
        cmat(j+k)=factor*1.d-3*d0_const*b_unit(k)
      ELSE IF (i .EQ. 4) THEN
        cmat(j+k)=factor*1.d-9*x_unit(k)*d_x
        !cmat(j+k)=factor*1.d-3*d0_const*x_unit(k)*d_x
        !cmat(j+k)=factor*1.d-3*d0_const*x_unit(k)*DSIN(u-u0)
      ELSE IF (i .EQ. 5) THEN
        cmat(j+k)=factor*1.d-3*d0_const*x_unit(k)*DSIN(3*u-u0)
      ELSE IF (i .EQ. 6) THEN
        cmat(j+k)=factor*1.d-3*d0_const*z_unit(k)*DSIN(u-u0)
      ELSE IF (i .EQ. 7) THEN
        cmat(j+k)=factor*1.d-3*d0_const*lambda*d_unit(k)*DCOS(u)
      ELSE IF (i .EQ. 8) THEN
        cmat(j+k)=factor*1.d-3*d0_const*lambda*d_unit(k)*DSIN(u)
      ELSE IF (i .EQ. 9) THEN
        cmat(j+k)=factor*1.d-3*d0_const*y_unit(k)*DCOS(u)
      ELSE IF (i .EQ. 10) THEN
        cmat(j+k)=factor*1.d-3*d0_const*y_unit(k)*DSIN(u)
      ELSE IF (i .EQ. 11) THEN
        cmat(j+k)=factor*1.d-3*d0_const*b_unit(k)*DCOS(u)
      ELSE IF (i .EQ. 12) THEN
        cmat(j+k)=factor*1.d-3*d0_const*b_unit(k)*DSIN(u)
      ELSE IF (i .EQ. 13) THEN
        cmat(j+k)=factor*1.d-3*d0_const*lambda*d_unit(k)*DCOS(2*u)
      ELSE IF (i .EQ. 14) THEN
        cmat(j+k)=factor*1.d-3*d0_const*lambda*d_unit(k)*DSIN(2*u)
      ELSE IF (i .EQ. 15) THEN
        cmat(j+k)=factor*1.d-3*d0_const*lambda*d_unit(k)*DCOS(4*u)
      ELSE IF (i .EQ. 16) THEN
        cmat(j+k)=factor*1.d-3*d0_const*lambda*d_unit(k)*DSIN(4*u)
      END IF
    END DO
  END DO
 

  !if (cprns .ne. trim(cprn)) then
  !  lfirst=.true.
  !  cprns=cprn
  !end if
  !lacc=.false.
  !if (lfirst .eq. .true.) then
  !  lfirst=.false.
  !  mjds=mjd
  !  sods=sod
  !  lacc=.true.
  !else
  !  if ((mjd-mjds)*86400+(sod-sods-60.0) .eq. 0.d0) then
  !     lacc=.true.
  !     mjds=mjd
  !     sods=sod
  !  end if
  !end if
  !if (lacc .eq. .true.) then
  !DO i=1, 3
  !  r_unit(i)=xsat(i)
  !  a_unit(i)=xsat(3+i)
  !  f(i)=f(i)+fp(i)*1.d3
  !  f(i)=f(i)/factor!-f_acr(i)
  !END DO
  !
  !CALL unit_vector(3,r_unit,r_unit,det)
  !CALL unit_vector(3,a_unit,a_unit,det)
  !CALL cross(r_unit,a_unit,c_unit)
  !CALL unit_vector(3,c_unit,c_unit,det)
  !CALL cross(c_unit,r_unit,a_unit)
  !CALL unit_vector(3,a_unit,a_unit,det)
  !f_acr(1)=dot(3,f,a_unit)
  !f_acr(2)=dot(3,f,c_unit)
  !f_acr(3)=dot(3,f,r_unit)
  !
  !f_xyz(1)=dot(3,f,x_unit)
  !f_xyz(2)=dot(3,f,y_unit)
  !f_xyz(3)=dot(3,f,z_unit)
  !
  !f_dyb(1)=dot(3,f,d_unit)
  !f_dyb(2)=dot(3,f,y_unit)
  !f_dyb(3)=dot(3,f,b_unit)
  !
  !WRITE(3000,'(I6,1X,3(F12.6,1X),A3,1X,21F12.6)') mjd,sod,beta*180.0/PI,u*180.0/PI,cprn,&
  !                   f*1.0d9,f_acr*1.0d9,f_dyb*1.0d9,f_xyz*1.0d9,DSQRT(f(1)**2+f(2)**2+f(3)**2)*1.0d9,factor
  
  !end if

  RETURN

END SUBROUTINE

