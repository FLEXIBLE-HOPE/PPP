!*
SUBROUTINE oi_solid_tides(PL,force_model,mjd,gmst,dc,ds)
!!
!! To compute the effect of solid Earth tides according to
!! IERS Conventions (2010), P82
!! Anelastic Earth model implemented
!! the total tidal contribution is computed, including the
!! time independent (permanent) contribution to the geopotential coefficient C20,
!! which is adequate for a "conventional tide free" model
!!
!*
USE orbit
IMPLICIT NONE

!*
! The arguments
!!--------------------------
TYPE(PLANET_INFO) :: PL
CHARACTER(LEN=*) :: force_model
REAL(RL) :: mjd,gmst
REAL(RL) :: dc(4,0:4),ds(4,0:4)

  !*
  ! The local variables
  !!--------------------------
  TYPE(SOLID_TIDE_FREQ) :: estf(0:2)

  INTEGER(IT) :: n,m,i

  REAL(RL) :: sinth,costh,theta

  REAL(RL) :: rek(2:3,0:3)
  REAL(RL) :: imk(2:3,0:3)
  REAL(RL) :: k2m(0:2)
  REAL(RL) :: P(3,0:3)
  REAL(RL) :: norm(3,0:3)
  REAL(RL) :: dcf(0:2),dsf(0:2)
  REAL(RL) :: arg(5),beta(6)
  REAL(RL) :: slat,lon,rho,fac

  LOGICAL(LG) :: lfirst
  DATA lfirst /.TRUE./

  DATA rek / &
   0.30190d0, 0.093d0, &
   0.29830d0, 0.093d0, &
   0.30102d0, 0.093d0, &
   0.00000d0, 0.094d0 /

  DATA imk / &
  -0.00000d0, 0.000d0, &
  -0.00144d0, 0.000d0, &
  -0.00130d0, 0.000d0, &
   0.00000d0, 0.000d0 /

  DATA k2m /-0.00089d0,-0.00080d0,-0.00057d0/

  SAVE lfirst,norm,rek,imk,k2m,estf

  !*
  ! Start the exectuable code
  !!--------------------------

  IF (lfirst .EQ. .TRUE.) THEN
    lfirst=.FALSE.
    norm=0.d0
    CALL normalize_coeff('norm',3,3,3,3,norm)
    CALL read_tide_freq('IERS2010',estf)
    DO n=2, 3
      DO m=0, n
        norm(n,m)=norm(n,m)/(2*n+1.d0)
      END DO
    END DO
  END IF

  !! Step 1
  DO i=1, PL.nplanet
    IF (INDEX(force_model,PL.name(i)(1:3)) .EQ. 0) CYCLE
    slat=PL.xe(3,i)/PL.dist2cb(i)
    lon =DATAN2(PL.xe(2,i),PL.xe(1,i))
    rho =PL.radius(PL.icb)/PL.dist2cb(i)
    fac =PL.gm(i)/PL.gm(PL.icb)
    P(2,0)=0.5d0*(3.d0*slat*slat-1.d0)
    P(2,1)=3.d0*DSQRT(1.d0-slat*slat)*slat
    P(2,2)=3.d0*(1.d0-slat*slat)
    P(3,0)=0.5d0*(5.d0*slat*slat-3.d0)*slat
    P(3,1)=1.5d0*DSQRT(1.d0-slat*slat)*(5.d0*slat*slat-1.d0)
    P(3,2)=3.d0*5.d0*(1.d0-slat*slat)*slat
    P(3,3)=3.d0*5.d0*DSQRT((1.d0-slat*slat)**3)

    ! IERS 2010, P82, Equ.6.6
    DO n=2, 3
      DO m=0, n
        dc(n,m)=dc(n,m)+norm(n,m)*fac*P(n,m)*(rek(n,m)*DCOS(m*lon)+imk(n,m)*DSIN(m*lon))*rho**(n+1)
        ds(n,m)=ds(n,m)+norm(n,m)*fac*P(n,m)*(rek(n,m)*DSIN(m*lon)-imk(n,m)*DCOS(m*lon))*rho**(n+1)
      END DO
    END DO

    ! IERS 2010, P83, Equ.6.7
    DO m=0, 2
      dc(4,m)=dc(4,m)+k2m(m)*norm(2,m)*fac*P(2,m)*DCOS(m*lon)*rho**3
      ds(4,m)=ds(4,m)+k2m(m)*norm(2,m)*fac*P(2,m)*DSIN(m*lon)*rho**3
    END DO
  END DO

  !! step 2
  dcf=0.d0
  dsf=0.d0

  CALL fund_arg_nutation(mjd,arg)
  CALL doodson_arg(gmst,arg,beta)

  DO m=0, 2
    DO i=1, estf(m).n
      theta=estf(m).etf(1,i)*beta(1)+estf(m).etf(2,i)*beta(2)+estf(m).etf(3,i)*beta(3)+ &
            estf(m).etf(4,i)*beta(4)+estf(m).etf(5,i)*beta(5)+estf(m).etf(6,i)*beta(6)
      sinth=DSIN(theta)
      costh=DCOS(theta)
      SELECT CASE(m)
        CASE(0)
          dcf(m)=dcf(m)+estf(m).amp(1,i)*costh-estf(m).amp(2,i)*sinth
          dsf(m)=0.d0
        CASE(1)
          dcf(m)=dcf(m)+estf(m).amp(1,i)*sinth+estf(m).amp(2,i)*costh
          dsf(m)=dsf(m)+estf(m).amp(1,i)*costh-estf(m).amp(2,i)*sinth
        CASE(2)
          dcf(m)=dcf(m)+estf(m).amp(1,i)*costh-estf(m).amp(2,i)*sinth
          dsf(m)=dsf(m)-estf(m).amp(1,i)*sinth-estf(m).amp(2,i)*costh
      END SELECT
    END DO
    dcf(m)=dcf(m)*1.d-12
    dsf(m)=dsf(m)*1.d-12
  END DO

  ! frequency-dependent computation makes dsf(0) non-zero,
  ! though the value is rather small,
  ! according to definition of spherical harmonics, it should be zero.
  !! need to be test
  !dsf(0)=0.d0

  !! step 3 used for 'zero tide' geopotential model to remove the permanent tide
  !! dcf(0)=dcf(0)-4.4228.d-8*(-0.31460)*0.29525d0

  DO m=0, 2
    dc(2,m)=dc(2,m)+dcf(m)
    ds(2,m)=ds(2,m)+dsf(m)
  END DO

  RETURN

END SUBROUTINE
