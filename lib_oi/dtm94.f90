!*
SUBROUTINE dtm94(mjd,gsat,gsun,f107,f107a,kp,den)
!!
!! COMPUTE THE ATMOSPHERE DENSITY AT THE SPECIFICAL SATELLITE POSITION
!! AND EPOCH ACCORDING TO DTM94 THERMOSPHERIC MODEL
!!
!! INPUT:
!!        MJD - MODIFIED JULIDAY DAY
!!        GSAT - GEODETIC COORDINATES (LATITUDE, LONGITUDE, ALTITUDE)
!!               OF THE SATELLITE (METERS)
!!        GSUN - GEODETIC COORDINATES (LATITUDE, LONGITUDE, ALTITUDE)
!!               OF THE SUN (METERS)
!!        KP - THE THREE-HOURLY GEOMAGNETIC INDEX TAKEN WITH A DELAY
!!             DEPENDING ON LATITUDE (3H AT THE POLE, 6H AT THE EQUATOR
!!             WITH A LINEAR INTERPOLATION)
!!        F107 - THE SOLAR RADIO FLUX AT 10.7 CM ON THE PREVIOUS DAY
!!        F107A - THE AVERAGE F107 OVER THREE SOLAR ROTATIONS (81 DAYS)
!!                BEFORE THE REQUIRED DAY
!!
!! OUTPUT:
!!        DEN - THERMOSPHERIC DENSITY (KG/M3)
!!
!! REFEREMCE:
!!      Berge C., Biancale R., Ill M., Barlier F., Improvement of the
!!        empirical thermoshpheric model DTM: DTM94-a comparative review
!!        of various temporal variations and prospects in space geodesy
!!        application, Journal of Geodesy,  Vol. 72 pp. 161-178 , 1998
!!
!!      Barlier F., Berger C., Falin JL, Kockarts G, Thuillier G,
!!        A thermospheric model based on satellite drag data, Annual
!!        Geophysicae, Vol. 34, pp. 9-24, 1978
!!
!!      Hedin AE, MSIS-86 thermospheric model. Journal of Geophysical
!!        Research, Vol. 92, pp. 4649-4662
!!
!!      Hedin AE et al, Empirical model of global thermospheric temperature
!!        and composition based on data from the Ogo-6 quadrupole mass
!!        spectrometer, Journal of Geophical Research, Vol 79, pp. 215-225, 1977
!!
!! AUTHOR: MAORONG GE, June-2003
!!
!*
USE const
IMPLICIT NONE

!*
! The arguments
!!-----------------------------
REAL(RL) :: mjd,gsat(1:*),gsun(1:*)
REAL(RL) :: kp,f107,f107a,den

  !*
  ! The local variables
  !!-----------------------------

  REAL(RL) :: t,d,g,t8,tz,pnm(12)
  REAL(RL) :: sigma,kusi,gamma,rho,fi

  INTEGER(IT) :: i,iyear,idoy

  LOGICAL(LG) :: lhohe

  ! model data
  REAL(RL) :: a(39,5),o2,mi(5),uamc,alpha(5)
  REAL(RL) :: tgrd,t120,h120,g120,boltman

  ! temperature
  DATA (a(i,5),i=1,39) &
       / 0.1000E+04, 0.9461E-02, 0.4267E-01, 0.1795E-02,-0.7990E-05, &
         0.3367E-02, 0.2263E-01, 0.3786E-01,-0.1923E-01,-0.9241E-02, &
        -0.2107E+03, 0.1032E-01, 0.2886E-01,-0.7631E+02,-0.1850E+00, &
        -0.2031E-01, 0.1447E-01,-0.3631E+01,-0.2894E-01,-0.1737E+03, &
        -0.1082E+00,-0.1999E-02, 0.3397E-02,-0.1613E-01,-0.9609E-02, &
        -0.1045E+00, 0.4575E-02, 0.4584E-02, 0.1776E-01,-0.4227E-02, &
        -0.3567E-02,-0.3606E-03, 0.1049E-01, 0.4571E-02,-0.2175E-04, &
         0.1511E-02, 0.2167E-02,-0.1142E-05, 0.6582E-04/
  ! H
  DATA (a(i,1),i=1,39) &
       / 0.1761E+06,-0.1337E+00, 0.0000E+00,-0.1246E-01, 0.0000E+00, &
        -0.1930E-01,-0.6000E-01,-0.2000E-01, 0.5878E-01, 0.0000E+00, &
         0.9227E+02, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.3301E+00, &
         0.1045E+00, 0.0000E+00,-0.1477E+02,-0.9065E-01,-0.7200E+02, &
         0.2094E+00, 0.2830E-01, 0.0000E+00, 0.8571E-01,-0.2475E-01, &
         0.3830E+00, 0.2941E-01, 0.0000E+00,-0.3974E-02, 0.4356E-01, &
         0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
         0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00/
  ! He
  DATA (a(i,2),i=1,39) &
       / 0.2791E+08, 0.1096E+00,-0.1908E+00,-0.2077E-03, 0.4835E-05, &
         0.2112E-02, 0.2212E-03,-0.1617E+00,-0.9222E-01,-0.8301E-02, &
         0.2135E+03, 0.2350E+00,-0.7905E-01, 0.1104E+03,-0.1268E+01, &
        -0.1851E+01, 0.6662E-01,-0.1870E+03,-0.4215E-01,-0.2167E+03, &
        -0.1278E+00,-0.6182E-02,-0.1745E-01,-0.4360E-01,-0.5285E-01, &
         0.3120E+00,-0.2137E-01,-0.2453E-01,-0.3673E-02,-0.8833E-01, &
         0.3399E-01, 0.5131E-02,-0.1474E-01, 0.8556E-02, 0.1957E-02, &
        -0.4503E-02,-0.1024E+00,-0.1687E-04,-0.1426E-02/
  ! O
  DATA (a(i,3),i=1,39) &
       / 0.8472E+11,-0.6645E-01,-0.9741E-01, 0.1228E-02, 0.4498E-05, &
         0.5358E-02, 0.2557E-02,-0.9822E-01, 0.1009E+00, 0.6261E-02, &
         0.1156E+02, 0.1760E+00,-0.7128E-01, 0.1064E+03, 0.3329E+00, &
        -0.1145E+00,-0.4242E-02,-0.6751E+00,-0.4171E-01, 0.1344E+03, &
        -0.6593E-01,-0.1934E-01,-0.8585E-02, 0.8635E-01, 0.8568E-01, &
         0.4549E-01,-0.3872E-01,-0.1287E-01,-0.9178E-01, 0.4262E-01, &
         0.4339E-01, 0.5689E-02, 0.6651E-02,-0.1017E-01, 0.2706E-02, &
        -0.4710E-02,-0.1464E-01,-0.6406E-05, 0.1518E-02/
  ! N2
  DATA (a(i,4),i=1,39) &
       / 0.3204E+12,-0.1402E+00, 0.5722E-01, 0.1126E-02,-0.2078E-05, &
         0.3424E-02,-0.1159E-01, 0.5516E-01,-0.1005E-01,-0.4390E-01, &
         0.1959E+03, 0.3248E-01, 0.5650E-01, 0.8820E+02, 0.2881E+00, &
        -0.3143E-01, 0.0000E+00,-0.2001E+03, 0.6070E-01, 0.5115E+02, &
        -0.4610E-01,-0.9700E-02, 0.3491E-02, 0.0000E+00, 0.0000E+00, &
        -0.7333E-01, 0.2116E-01, 0.7033E-02, 0.0000E+00, 0.0000E+00, &
        -0.6920E-02, 0.0000E+00,-0.7290E-02, 0.0000E+00,-0.2920E-02, &
         0.2705E-02,-0.5994E-02, 0.0000E+00, 0.8213E-02/

  ! O2
  DATA o2 / 0.4475E+11 /

  ! molecular mass
  DATA mi /2.E0, 4.E0, 16.E0, 28.E0, 32.E0/
  !DATA mi /1.0079E0, 4.0026E0, 15.999E0, 28.0134E0, 31.998E0/

  ! Atomic mass unit (kg) 1.660538 921(73)×10-27 kg（CODATA2010）
  DATA uamc /1.6605402E-27/

  ! thermal diffusion coefficient for H and He
  DATA alpha /-0.38E0,-0.38E0,0.E0,0.E0,0.E0/

  ! unit: K/Km K KM m/s/s ,  N*m*K  Km
  DATA tgrd,t120,h120,g120 /14.348E0,380.0E0,120.0E0,9.446626E0/
  DATA boltman /1.380310E-23/

  !*
  ! Start the exectuable code
  !!-----------------------------

  ! local sideral time in hours
  t=DMOD((gsat(2)-gsun(2))/PI*12.d0+12.d0,24.d0)
  IF (t .LT. 0) t=t+24.0

  ! day of year in day
  CALL mjd2doy(INT(mjd),iyear,idoy)
  d=idoy+(mjd-INT(mjd))

  CALL dtm94_legendre(gsat(1),pnm)

  ! G function for temperature
  CALL dtm94_gfunction(.TRUE.,.FALSE.,a(1,5),pnm,f107,f107a,kp,d,t,g)
  ! thermopause temperature in K
  t8=a(1,5)*(1.0+g)

  ! relative vertical temperature gradient T' in Km-1
  sigma=tgrd/(t8-t120)

  ! geopotential altitude in Km
  kusi=(6356.770+h120)/(6356.770+gsat(3)*1.d-3)*(gsat(3)*1.d-3-h120)
  ! temperature at the altitude z above the standard ellipsoid in K
  tz=t8-(t8-t120)*EXP(-sigma*kusi)

  rho=0.0
  DO i=1, 5
    lhohe=.TRUE.
    IF (i .GT. 2) lhohe=.FALSE.
    ! kg*m/s/s*km/(kg*m*m/s/s)/K=1000
    gamma=mi(i)*uamc*g120/sigma/boltman/t8*1.d3
    fi=(t120/tz)**(1+alpha(i)+gamma)*EXP(-sigma*gamma*kusi)

    IF (i .NE. 5) THEN
      CALL dtm94_gfunction(.FALSE.,lhohe,a(1,i),pnm,f107,f107a,kp,d,t,g)
      g=a(1,i)*EXP(g)
    ELSE
      g=o2
    END IF
    ! kg
    rho=rho+mi(i)*g*fi
  END DO

  ! kg/cm**3   ==> kg/m**3
  den=DBLE(rho*uamc*1.d6)

  RETURN

END SUBROUTINE


!*
SUBROUTINE dtm94_legendre(lat,pnm)
!!
!!
USE par
IMPLICIT NONE

!*
! The arguments
!!-----------------------------
REAL(RL) :: lat
REAL(RL) :: pnm(1:*)

  !*
  ! The local variables
  !!-----------------------------
  REAL(RL) :: sinlat,sinlat2,sinlat4,coslat,coslat2

  !*
  ! Start the exectuable code
  !!-----------------------------

  sinlat=DSIN(lat)
  sinlat2=sinlat*sinlat
  sinlat4=sinlat2*sinlat2
  coslat=DCOS(lat)
  coslat2=coslat*coslat

  ! p10,p20,p30,p40,p50
  pnm( 1)=sinlat
  pnm( 2)=0.5*(3.0*sinlat2-1.0)
  pnm( 3)=0.5*(5.0*sinlat2-3.0)*sinlat
  pnm( 4)=0.125*(35.0*sinlat4-30.0*sinlat2+3.0)
  pnm( 5)=0.125*(63.0*sinlat4-70.0*sinlat2+15.0)*sinlat

  ! p11,p21,p31,p51
  pnm( 6)=coslat
  pnm( 7)=3.0*sinlat*coslat
  pnm( 8)=0.5*(15.0*sinlat2-3.0)*coslat
  pnm( 9)=0.125*(315.0*sinlat4-210.0*sinlat2+15.0)*coslat

  ! p22,p32,p33
  pnm(10)=3.0*coslat2
  pnm(11)=15.0*sinlat*coslat2
  pnm(12)=15.0*coslat2*coslat

  RETURN

END SUBROUTINE


!*
SUBROUTINE dtm94_gfunction(ltem,lhohe,a,pnm,f107,f107a,kp,d,t,g)
!!
!!
!*
USE const
IMPLICIT NONE


!*
! The arguments
!!-----------------------------
LOGICAL(LG) :: ltem,lhohe
REAL(RL) :: a(1:*),pnm(1:*),f107,f107a,kp
REAL(RL) :: d,t,g

  !*
  ! The local variables
  !!-----------------------------
  INTEGER(IT) :: i
  REAL(RL) :: omg,omega
  REAL(RL) :: fact,solat

  !*
  ! Start the exectuable code
  !!-----------------------------
  omg=PI/12.0
  omega=2.0*PI/365.0

  solat=a(4)*(f107-f107a)+a(5)*(f107-f107a)**2+a(6)*(f107a-150.0)+a(38)*(f107a-150.0)**2
  IF (lhohe .EQ. .TRUE.) THEN
    fact=1.d0
  ELSE
    fact=1+solat
  END IF

  ! non periodic terms
  g=a(2)*pnm(2)+a(3)*pnm(4)+a(37)*pnm(1) &
    ! solar activity
    +solat &
    ! magnetic activity
    +(a(7)+a(8)*pnm(2))*kp &
    ! periodic terms
    +fact*( &
       ! symmetrical annual
        (a(9)+a(10)*pnm(2))*DCOS(omega*(d-a(11))) &
       ! symmetrical semi-annual
       +(a(12)+a(13)*pnm(2))*DCOS(2*omega*(d-a(14))) &
       ! asymmetrical annual (seasonal)
       +(a(15)*pnm(1)+a(16)*pnm(3)+a(17)*pnm(5))*DCOS(omega*(d-a(18))) &
       ! asymmetrical semi-annual
       +a(19)*pnm(1)*DCOS(2*omega*(d-a(20))) &
       ! diurnal
       +(a(21)*pnm(6)+a(22)*pnm(8)+a(23)*pnm(9)+(a(24)*pnm(6)+a(25)*pnm(7))*DCOS(omega*(d-a(18))))*DCOS(omg*t) &
       +(a(26)*pnm(6)+a(27)*pnm(8)+a(28)*pnm(9)+(a(29)*pnm(6)+a(30)*pnm(7))*DCOS(omega*(d-a(18))))*DSIN(omg*t) &
       ! semidiurnal
       +(a(31)*pnm(10)+a(32)*pnm(11)*DCOS(omega*(d-a(18))))*DCOS(2*omg*t) &
       +(a(33)*pnm(10)+a(34)*pnm(11)*DCOS(omega*(d-a(18))))*DSIN(2*omg*t) &
       ! terdiurnal
       +a(35)*pnm(12)*DCOS(3*omg*t)+a(36)*pnm(12)*DSIN(3*omg*t))


  ! for the temperature
  IF (ltem .EQ. .TRUE.) THEN
    g=g+a(39)*EXP(kp)
  ! for the neutral constituents
  ELSE
    g=g+a(39)*kp*kp
  END IF

  RETURN

END SUBROUTINE
