!*
SUBROUTINE q2rotmat(sc,q,rotmat)
!!
!! CALCULATE THE ROTATIOM MATRIX FROM QUANTERNION
!! SPACECRAFT-FIXED TO INTERTIAL SYSTEM
!*
USE par
IMPLICIT NONE

!*
! The arguments
!!-----------------------------
CHARACTER(LEN=*) :: sc
REAL(RL) :: q(4), rotmat(3,3)

  !*
  ! The local variables
  !!---------------------------
  REAL(RL) :: q11,q22,q33,q44,q12,q13,q14,q23,q24,q34

  !*
  ! Start the exectuable code
  !!---------------------------

  q11=q(1)*q(1)
  q22=q(2)*q(2)
  q33=q(3)*q(3)
  q44=q(4)*q(4)
  q12=q(1)*q(2)
  q13=q(1)*q(3)
  q14=q(1)*q(4)
  q23=q(2)*q(3)
  q24=q(2)*q(4)
  q34=q(3)*q(4)

  ! nominal definitation of Quaternion-derived rotation matrix
  IF (sc(1:5) .EQ. 'JASON') THEN
    rotmat(1,1)=q11+q22-q33-q44
    rotmat(1,2)=2.d0*(q23-q14)
    rotmat(1,3)=2.d0*(q24+q13)
    rotmat(2,1)=2.d0*(q23+q14)
    rotmat(2,2)=-q44-q22+q11+q33
    rotmat(2,3)=2.d0*(q34-q12)
    rotmat(3,1)=2.d0*(q24-q13)
    rotmat(3,2)=2.d0*(q34+q12)
    rotmat(3,3)=q44-q33+q11-q22
  ELSE IF (sc(1:5) .EQ. 'SWARM') THEN
    rotmat(1,1)=1.d0-2.d0*q22-2.d0*q33
    rotmat(1,2)=2.d0*(q12+q34)
    rotmat(1,3)=2.d0*(q13-q24)
    rotmat(2,1)=2.d0*(q12-q34)
    rotmat(2,2)=1.d0-2.d0*q11-2.d0*q33
    rotmat(2,3)=2.d0*(q23+q14)
    rotmat(3,1)=2.d0*(q13+q24)
    rotmat(3,2)=2.d0*(q14-q23)
    rotmat(3,3)=1.d0-2.d0*q11-2.d0*q22
  ELSE
    rotmat(1,1)=q44+q11-q22-q33
    rotmat(1,2)=2.d0*(q12-q34)
    rotmat(1,3)=2.d0*(q13+q24)
    rotmat(2,1)=2.d0*(q12+q34)
    rotmat(2,2)=q44+q22-q11-q33
    rotmat(2,3)=2.d0*(q23-q14)
    rotmat(3,1)=2.d0*(q13-q24)
    rotmat(3,2)=2.d0*(q23+q14)
    rotmat(3,3)=q44+q33-q11-q22
  END IF

  RETURN

END SUBROUTINE
