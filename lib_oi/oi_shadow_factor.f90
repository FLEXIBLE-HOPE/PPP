!*
SUBROUTINE oi_shadow_factor(rad_sun,rad_earth,rad_moon,xsun,xsat,xmoon,lambda)
!!
!*
USE const
IMPLICIT NONE

!*
! The arguments
!!-----------------------------
REAL(RL) :: rad_sun,rad_earth,rad_moon
REAL(RL) :: xsun(1:*),xsat(1:*),xmoon(1:*),lambda

  !*
  ! The local variables
  !!-----------------------------
  INTEGER(IT) :: i
  REAL(RL) :: dt,sep,dsun,dearth,dmoon
  REAL(RL) :: lsat2sun,lsat2earth,lsun2earth,lmoon2sun,lsat2moon
  REAL(RL) :: u_xsat2earth(3),u_xsat2sun(3),xmoon2sun(3),u_xsat2moon(3)

  !*
  ! The functions called
  !!-----------------------------
  REAL(RL) :: dot


  !*
  ! Start the exectuable codes
  !!-----------------------------

  ! Get sight angle (arc) between the sun and the earth seen from satellite
  lsun2earth=0.d0
  lsat2sun=0.d0
  lsat2earth=0.d0
  lmoon2sun=0.d0
  lsat2moon=0.d0
  lambda=1.0d0
  
  DO i=1, 3
    lsun2earth=lsun2earth+xsun(i)*xsun(i)
    lsat2sun=lsat2sun+(xsun(i)-xsat(i))**2
  END DO
  lsun2earth=DSQRT(lsun2earth)
  lsat2sun=DSQRT(lsat2sun)
  
  dt=(lsat2sun-lsun2earth)/VEL_LIGHT
  DO i=1, 3
    u_xsat2sun(i)=xsun(i)-xsat(i)-xsat(i+3)*dt
    u_xsat2earth(i)=-xsat(i)-xsat(i+3)*dt
  END DO

  CALL unit_vector(3,u_xsat2sun,u_xsat2sun,lsat2sun)

  IF(lsat2sun .LE. lsun2earth) GOTO 200
  
  CALL unit_vector(3,u_xsat2earth,u_xsat2earth,lsat2earth)

  sep=dot(3,u_xsat2sun,u_xsat2earth)
  sep=DACOS(sep)
  
  dsun=rad_sun/lsat2sun
  dearth=rad_earth/lsat2earth

  ! sight angle (arc) between for sun and earth
  CALL shadow_factor(dsun,dearth,sep,lambda)
  
  ! If no earth eclipse, check the moon
200 CONTINUE
  IF (lambda .LT. 1.d0) THEN
    RETURN
  ELSE
    DO i=1,3
      xmoon2sun(i)=xmoon(i)-xsun(i)
    END DO
            
    lmoon2sun=DSQRT(xmoon2sun(1)**2+xmoon2sun(2)**2+xmoon2sun(3)**2)
    dt=(lsat2sun-lmoon2sun)/VEL_LIGHT
    DO i=1,3
      u_xsat2moon(i)=xmoon(i)-xsat(i)-xsat(i+3)*dt
    END DO

    IF(lsat2sun .LE. lmoon2sun) RETURN
      
    CALL unit_vector(3,u_xsat2moon,u_xsat2moon,lsat2moon)
    dsun=rad_sun/lsat2sun
    dmoon=rad_moon/lsat2moon
      
    sep=dot(3,u_xsat2sun,u_xsat2moon)
    sep=DACOS(sep)

    CALL shadow_factor(dsun,dmoon,sep,lambda)
  END IF

  RETURN
      
END SUBROUTINE


SUBROUTINE shadow_factor(rs,rp,sep,lambda)
!!
!*
USE const
IMPLICIT NONE
 
!*
! The arguments
!!-----------------------------
REAL(RL) :: rs,rp,sep,lambda

  !*
  ! The local variables
  !!-------------------------------
  REAL(RL) :: r1,r2,phi,ari
  REAL(RL)::  hgt,thet,area1,area2,area3

  !*
  ! The function called
  !!-------------------------------
  REAL(RL) :: dot

  !*
  ! Start the exectuable codes
  !!-------------------------------

  IF (rs+rp .LE. sep) THEN
    ! no eclipse
    RETURN
  ELSE IF (rp-rs .GT. sep) THEN
    ! full eclipse
    lambda=0.d0
    RETURN
  ELSEIF(sep .LE. rs-rp) THEN
    ! partial eclipse, do the calculations
    lambda=(rs**2-rp**2)/rs**2
  ELSE
    ! set r1=smaller disc, r2=larger
    r1=MIN(rp,rs)
    r2=MAX(rp,rs)
       
    phi=DACOS((r1*r1+sep*sep-r2*r2)/(2.0d0*r1*sep))
    IF (phi .LT. 0.d0) phi=PI+phi
    IF (r2/r1 .GT. 5.0d0) THEN
      hgt=DSQRT(r1**2-(sep-r2)**2)
      area2=hgt*(sep-r2)
      area3=0.0d0
    ELSE
      hgt=r1*DSIN(phi)
      thet = DASIN(hgt/r2)
      area2=sep*hgt
      area3=thet*r2**2
    END IF
     
    area1= (PI-phi)*r1**2
    ari=area1+area2-area3
    area1=PI*rs**2
    IF (rs .GT. rp) THEN
      area2=PI*rp**2
      lambda=(area1+ari-area2)/area1
    ELSE
      lambda = ari/area1
    END IF
  END IF
      
  RETURN
    
END SUBROUTINE
