!*
SUBROUTINE brd2xyz(cprn,neph,eph,wk,sow,xsat,clk)
!!
!! To be carefully, different navigation message has been broadcasted by GNSS system
!! Type 1: QZSS (LNAV, L1C/A), IRNSS
!! Type 2: QZSS (CNAV2, L1C; CNAV, L2C/L5)
!!
!*
USE brdeph
IMPLICIT NONE

!*
! The arguments
!!---------------------
ChARACTER(LEN_PRN) :: cprn
TYPE(GPS_BRDEPH) :: eph(1:*)
INTEGER(IT) :: neph,wk
REAL(RL) :: sow,clk,xsat(1:*)

  !*
  ! The local variables
  !!-------------------------

  INTEGER(IT) :: i,ieph
  REAL(RL) :: dt,dtc,dtmin,pos(3)
  REAL(RL) :: GMEARTH,EARTH_ROTATE
  REAL(RL) :: a,xn,xm,ex,e,v0,vs,vc,phi,ccc,sss,du,dr,di,r,u,xi,xx,yy
  REAL(RL) :: xnode,term,xpdot,ypdot,asc,xinc,xp,yp,asctrm,v,clkr

  !*
  ! Start of executable code
  !!------------------------

  SELECT CASE(cprn(1:1))
    CASE('G')
      GMEARTH = 3.986005D14
      EARTH_ROTATE = 7.2921151467D-5
    CASE('E')
      GMEARTH = 3.986004418D14
      EARTH_ROTATE = 7.2921151467D-5
    CASE('C')
      GMEARTH = 3.986004418D14
      EARTH_ROTATE = 7.2921150D-5
    CASE('J')
      GMEARTH = 3.986005D14
      EARTH_ROTATE = 7.2921150D-5
    CASE('I')
      GMEARTH = 3.986005D14
      EARTH_ROTATE = 7.29211514670D-5
    CASE DEFAULT
      GMEARTH = GME
      EARTH_ROTATE = E_ROTATE
  END SELECT

  ieph=0
  dtmin=12.d0
  DO i=1, neph
    IF (eph(i).cprn .EQ. cprn) THEN
      dt=(wk-eph(i).week)*168.d0+(sow-eph(i).toe)/3600.d0
      IF(ABS(dt) .LE. dtmin) THEN
        dtmin=ABS(dt)
        ieph=i
      END IF
    END IF
  END DO
  IF(ieph .EQ. 0) RETURN

  dt=(wk-eph(ieph).week)*604800.d0+sow-eph(ieph).toe
  dtc=(wk*7+44244-eph(ieph).mjd)*86400.d0+sow-eph(ieph).sod
  clk=eph(ieph).a0+(eph(ieph).a1+eph(ieph).a2*dtc)*dtc

  a=eph(ieph).roota**2
  xn=DSQRT(GMEARTH/a/a/a)
  xn=xn+eph(ieph).dn

  xm=eph(ieph).m0+xn*dt
  ex=xm
  e=eph(ieph).e
  DO i=1, 12
    ex=xm+e*DSIN(ex)
  END DO

  ! This has been corrected in gpsmodel
  !clkr=-2.d0*DSQRT(GME)/VEL_LIGHT/VEL_LIGHT*eph(ieph).e*eph(ieph).roota*DSIN(ex)
  !clk=clk+clkr

  v0=1.d0-e*DCOS(ex)
  vs=DSQRT(1.d0-e*e)*DSIN(ex)/v0
  vc=(DCOS(ex)-e)/v0
  !v=DABS(DASIN(vs))
  !IF (vc .GE. 0.d0) THEN
  !  IF (vs .LT. 0.d0) v=2.d0*PI-v
  !ELSE
  !  IF(vs .LE. 0.d0) THEN
  !    v=PI+v
  !  ELSE
  !    v=PI-v
  !  END IF
  !END IF
  v=DATAN2(vs,vc)

  phi=v+eph(ieph).omega

  ccc=DCOS(2.d0*phi)
  sss=DSIN(2.d0*phi)
  du=eph(ieph).cuc*ccc+eph(ieph).cus*sss
  dr=eph(ieph).crc*ccc+eph(ieph).crs*sss
  di=eph(ieph).cic*ccc+eph(ieph).cis*sss
  r=a*(1.d0-e*DCOS(ex))+dr
  u=phi+du

  xi=eph(ieph).i0+di+eph(ieph).idot*dt
  xx=r*DCOS(u)
  yy=r*DSIN(u)

  IF (cprn.EQ.'C01' .OR. cprn.EQ.'C02' .OR. cprn.EQ.'C03' .OR. &
      cprn.EQ.'C04' .OR. cprn.EQ.'C05' .OR. cprn.EQ.'C17' .OR. &
      cprn.EQ.'C59') THEN
    xnode=eph(ieph).omega0+eph(ieph).omegadot*dt
  ELSE
    xnode=eph(ieph).omega0+(eph(ieph).omegadot-EARTH_ROTATE)*dt
  END IF
  !! if the time system is different, it is better to compute the orbit
  !! in their own system, because TOE will be different if you transform
  !! it to other system
  xnode=xnode-EARTH_ROTATE*eph(ieph).toe
  xsat(1)=xx*DCOS(xnode)-yy*DCOS(xi)*DSIN(xnode)
  xsat(2)=xx*DSIN(xnode)+yy*DCOS(xi)*DCOS(xnode)
  xsat(3)=yy*DSIN(xi)

  ! rx(-5.0*rad) for GEO
  IF (cprn.EQ.'C01' .OR. cprn.EQ.'C02' .OR. cprn.EQ.'C03' .OR. &
      cprn.EQ.'C04' .OR. cprn.EQ.'C05' .OR. cprn.EQ.'C17' .OR. &
      cprn.EQ.'C59') THEN
    pos(1)=xsat(1)
    pos(2)=DCOS(-5.0*DEG2RAD)*xsat(2)+DSIN(-5.0*DEG2RAD)*xsat(3)
    pos(3)=DSIN(5.0*DEG2RAD)*xsat(2)+DCOS(-5.0*DEG2RAD)*xsat(3)

    ! rz(wearth*dt)
    xsat(1)=pos(1)*DCOS(EARTH_ROTATE*dt)+pos(2)*DSIN(EARTH_ROTATE*dt)
    xsat(2)=-1.d0*pos(1)*DSIN(EARTH_ROTATE*dt)+pos(2)*DCOS(EARTH_ROTATE*dt)
    xsat(3)=pos(3)
  END IF

  term=(xn*a)/DSQRT(1.d0-e*e)
  xpdot=-DSIN(u)*term
  ypdot=(e+DCOS(u))*term
  asc=xnode
  xinc=xi
  xp=xx
  yp=yy

  ! The velocity is not right for GEO
  asctrm=(eph(ieph).omegadot-EARTH_ROTATE)
  xsat(4)=xpdot*DCOS(asc)-ypdot*DCOS(xinc)*DSIN(asc)-xp*DSIN(asc)*asctrm-yp*DCOS(xinc)*DCOS(asc)*asctrm
  xsat(5)=xpdot*DSIN(asc)+ypdot*DCOS(xinc)*DCOS(asc)+xp*DCOS(asc)*asctrm-yp*DCOS(xinc)*DSIN(asc)*asctrm
  xsat(6)=ypdot*DSIN(xinc)


  RETURN

END SUBROUTINE
