!*
SUBROUTINE phase_windup(lsit,lfirst,rot_f2j,rot_l2f,xbf,ybf,zbf,xrec2sat,dphi0,dphi)
!!
!! purpose  : phase wind-up correction ( the receiver and satellite antenna orientation dependent
!!            phase corrections). See Wu J.T., et al., Manuscripta Geogetica (1993) 18, pp91-98
!!
!! parameter:
!!           lfirst -- first call for this satellite-station pair
!!           rot_f2j -- rotation matrix from earth-fixed to inertial (J2000)
!!           rot_l2f -- rotation matrix from station system (enu right-hand system) to earth-fixed
!!           xbf,ybf,zbf -- unit vectors of spacecraft-fixed system (rotation matrix from body-fixed 
!!                          to inertial)
!!           xrec2sat -- vector from rec. to satellite  for k-direction
!!           dphi0 -- initial dphi
!!           dphi  --  phase correction
!!
!! author   : Maorong Ge
!!
!! modified : 5/26/2007 -- SIGN used for dphi
!!            5/26/2007 -- add leo antenna system
!!
!*
USE par
USE const
IMPLICIT NONE

!*
! The arguments
!!---------------------
LOGICAL(LG) :: lsit,lfirst
REAL(RL) :: rot_f2j(3,3),rot_l2f(3,3),xbf(3),ybf(3),zbf(3),xrec2sat(3)
REAL(RL) :: dphi0,dphi

  !*
  ! The local variables
  !!----------------------
  INTEGER(IT) :: j,n
  REAL(RL) :: x_l(3),y_l(3),x_f(3),y_f(3),x_j(3),y_j(3),dummy(3),rlength,kusi
  REAL(RL) :: k(3),d_r(3),d_s(3)

  !! enu (right-hand system) for sit
  DATA y_l/-1.0,0.0,0.0/

  !*
  ! The function called
  !!----------------------
  REAL(RL) :: dot

  !*
  ! Start the exectuable code
  !!--------------------------

  IF (lsit .EQ. .FALSE.) THEN
    x_l(1)= 0.d0
    x_l(2)=-1.d0
    x_l(3)= 0.d0
  ELSE
    x_l(1)= 0.d0
    x_l(2)= 1.d0
    x_l(3)= 0.d0
  END IF

  !! unit vector of the signal transmitting direction
  CALL unit_vector(3,xrec2sat,k,rlength)
  DO j=1,3
    k(j)=-k(j)
  END DO

  !
  !! equivalent antenna dipole for both receiver and satellite antenna
  !!    D = x_j - k (k . x_j) - k x y_j
  !! k is the unit vector of signal transmitting dirrection
  !
  !! transfer the unit vector of the antenna dipole unit in local/antenna system to
  !! inertial system.
  !! For ground receiver local   => earth fixed      => inertial
  !! For LEO    receiver antenna => spacecraft fixed => inertial
  CALL matmpy(rot_l2f,x_l,x_f,3,3,1)
  CALL matmpy(rot_l2f,y_l,y_f,3,3,1)
  CALL matmpy(rot_f2j,x_f,x_j,3,3,1)
  CALL matmpy(rot_f2j,y_f,y_j,3,3,1)

  !! D = x_j - k(k . x_j) + k x y_j
  rlength=dot(3,k,x_j)
  CALL cross(k,y_j,dummy)
  DO j=1,3
    d_r(j)=x_j(j)-rlength*k(j)+dummy(j)
  END DO
  CALL unit_vector(3,d_r,d_r,rlength)

  !! the same for the satellite antenna
  rlength=dot(3,k,xbf)
  CALL cross(k,ybf,dummy)
  DO j=1,3
    d_s(j)=xbf(j)-rlength*k(j)-dummy(j)
  END DO
  CALL unit_vector(3,d_s,d_s,rlength)

  !! kusi  k . (D_s x D_r)
  CALL cross(d_s,d_r,dummy)
  kusi=dot(3,k,dummy)
  dphi=dot(3,d_s,d_r)
  IF (DABS(dphi) .GE. 1.d0) dphi=DSIGN(1.d0,dphi)
  dphi=DSIGN(1.d0,kusi)*DACOS(dphi)/(2*PI)

  IF (lfirst) THEN
    lfirst=.FALSE.
    n=0
  ELSE
    n=NINT(dphi0-dphi)
  END IF
  dphi=n+dphi

  !! save for the next epoch
  dphi0=dphi

  RETURN

END SUBROUTINE
