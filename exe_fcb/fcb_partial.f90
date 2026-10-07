!! purpose  : observation equaiton of a GPS observation
!!
!! parameter: npar -- number of parameters of the station
!!            ur2s -- partial of range wrt satellite position.
!!            drate -- range rate 
!!            parname -- name of the local parameters
!!            ltog -- index of local parameter in global parameter table
!!            xrot    -- rotation matrix for earthfixed system to inertial system
!!            trpart  -- partial of range wrt ztd and htg
!!            amat    -- observation equaiton according to local parameter table  
!

SUBROUTINE fcb_partial(npar,ur2s,drate,parname,ltog,xrot,trpart,grdpart,amat)
!!
!*
USE par
IMPLICIT NONE

!*
! The arguments
!!-------------------------
INTEGER(IT) :: npar, ltog(MAXPARSIT)
REAL(RL) :: ur2s(3),drate,xrot(3,3),trpart,grdpart(2),amat(MAXPARSIT)
CHARACTER(LEN_PARNAME) :: parname(MAXPARSIT)

  !*
  ! The local variables
  !!-----------------------
  INTEGER(IT) :: i, ipar

  !*
  ! Start the exectuable code
  !!----------------------------

  DO ipar=1, MAXPARSIT
    amat(ipar)=0.d0
  END DO

  ipar=0
  DO WHILE(ipar .LT. npar)
    ipar=ipar+1
    IF (ltog(ipar) .EQ. 0) CYCLE

    !! station coordinates
    IF (parname(ipar)(1:5) .EQ. 'STAPX') THEN
      DO i=1,3
        amat(ipar+i-1) =-xrot(1,i)*ur2s(1)-xrot(2,i)*ur2s(2)-xrot(3,i)*ur2s(3)
      END DO
      ipar=ipar+2

    !! atmospheric delay
    ELSE IF (parname(ipar)(1:3) .EQ. 'ZTD') THEN
      amat(ipar)=trpart


    ELSE IF (parname(ipar)(2:4) .EQ. 'GRD') THEN
      amat(ipar)=grdpart(1)
      amat(ipar+1)=grdpart(2)
      ipar=ipar+1

    ELSE IF (parname(ipar)(1:3) .EQ. 'ION') THEN
      amat(ipar)=1.d0

    !! receiver clock epoch
    ELSE IF (parname(ipar)(1:6) .EQ. 'RECCLK') THEN
      amat(ipar)=1.d0 !-drate

    !! IFHB in receiver side
    ELSE IF (parname(ipar)(1:6) .EQ. 'RECDCB') THEN
      amat(ipar)=1.d0

    !! ambiguity parameters
    ELSE IF (parname(ipar)(1:3) .EQ. 'AMB') THEN
      amat(ipar)=1.d0
    END IF
  END DO

  RETURN

END SUBROUTINE
