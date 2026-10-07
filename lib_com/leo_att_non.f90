!*
SUBROUTINE leo_att_non(snam,xsat,rmat)
!!
!!  PURPOSE: COMPUTE THE ROTATION MATRIX FROM NOMINAL SPACECRAFT-FIXED SYSTEM
!!           TO INERTIAL SYSTEM
!!
!! NOTE :
!!      CHAMP   z -- to earth center  (U)
!!              x -- aligned with the long side of the sc towards the boom (E)
!!                   in flight direction
!!              y -- forming a right-handed-system with x and z (N)
!!
!!      GRACE  (SRF to INERTIAL SYSTEM)
!!              x -- from origin to KB/KBa phase center (E)
!!                   GRACE_A  -v,  GRACE_B +v
!!              z -- normal to x and to the plane of the main equipment platform (in nadir dire) (U)
!!              y -- forming a right-handed-system with x and z (N)
!!
!!      COSMIC  z -- nadir = -r/|r| (U)
!!              y -- (r x v) / | r x v |
!!              x -- y x z
!!      JASON   z -- nadir = -r/|r| (U)
!!              x -- v/|v|
!!              y -- z * x
!!
!*
USE const
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!-----------------------------
CHARACTER(LEN=*) :: snam
REAL(RL) :: xsat(1:*),rmat(3,3)

  !*
  ! The local variables
  !!-----------------------------

  INTEGER(IT) :: i,kmod

  REAL(RL) :: ux(3),uv(3),un(3)
  REAL(RL) :: rotm(3,3),lent

  !*
  ! Start the exectuable code
  !!----------------------------

  kmod=0

  SELECT CASE(TRIM(snam))
    CASE('GRACEB')
      kmod=-1
    CASE DEFAULT
      kmod=1
  END SELECT

  IF (kmod .EQ. 0) THEN
    WRITE(ERROR_UNIT,'(A)') '***ERROR(leo_att_non): unknown LEO '//TRIM(snam)
    CALL exit(1)
  END IF

  DO i=1, 3
    ux(i)=-xsat(i)
    uv(i)=xsat(3+i)
    IF (kmod .EQ. -1) uv(i)=-uv(i)
  END DO

  CALL unit_vector(3,ux,ux,lent)
  CALL unit_vector(3,uv,uv,lent)
  CALL cross(ux,uv,un)
  CALL unit_vector(3,un,un,lent)

  SELECT CASE(TRIM(snam))
    CASE('COSMIC','JASON','TECHDEMSAT1')
      CALL cross(un,ux,uv)
      CALL unit_vector(3,uv,uv,lent)
    CASE DEFAULT
      CALL cross(uv,un,ux)
      CALL unit_vector(3,ux,ux,lent)
  END SELECT

  DO i=1,3
    rmat(i,1)=uv(i)
    rmat(i,2)=un(i)
    rmat(i,3)=ux(i)
  END DO

  IF (TRIM(snam) .EQ. 'GAOFEN3A') THEN
    rotm=0.d0
    rotm(1,1)=1.d0
    !! left side viewing
    rotm(2,2)=DCOS(31.5*DEG2RAD)
    rotm(2,3)=DSIN(31.5*DEG2RAD)
    rotm(3,2)=-rotm(2,3)
    rotm(3,3)=rotm(2,2)
    !! right side viewing: the mian mode
    ! rotm(2,2)=DCOS(-31.5*DEG2RAD)
    ! rotm(2,3)=-DSIN(-31.5*DEG2RAD)
    ! rotm(3,2)=DSIN(-31.5*RAD2DEG)
    ! rotm(3,3)=DCOS(-31.5*DEG2RAD)
    CALL matmpy(rmat,rotm,rmat,3,3,3)
    rotm=0.d0
    rotm(2,2)=1.d0
    rotm(1,1)=DCOS(5.d0*DEG2RAD)
    rotm(1,3)=DSIN(5.d0*DEG2RAD)
    rotm(3,1)=-rotm(1,3)
    rotm(3,3)=rotm(1,1)
    CALL matmpy(rmat,rotm,rmat,3,3,3)
    rotm=0.d0
    rotm(3,3)=1.d0
    rotm(1,1)=DCOS(-7.d0*DEG2RAD)
    rotm(1,2)=DSIN(-7.d0*DEG2RAD)
    rotm(2,1)=-rotm(1,2)
    rotm(2,2)=rotm(1,1)
    CALL matmpy(rmat,rotm,rmat,3,3,3)
  END IF

  RETURN

END SUBROUTINE
