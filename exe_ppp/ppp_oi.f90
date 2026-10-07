!*
SUBROUTINE ppp_oi(mjd,sod,dintv,SAT,xsat)
!!
!*
USE orbit
USE satellite
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!------------------------
INTEGER(IT) :: mjd
REAL(RL) :: sod,dintv,xsat(6)
TYPE(SATE) :: SAT
  
  !*
  ! The local variables
  !!---------------------------
  TYPE(ORBCFG) :: CKF
  

  INTEGER(IT) :: i,j,k,kk,nsign
  INTEGER(IT) :: isat,idir,ii

  REAL(RL) :: sod0,sod1
  REAL(RL) :: tb,te,t,tdir,hh
  REAL(RL) :: x(MAXEQUS),acc(MAXEQUS)
  REAL(RL) :: amat(3,3),bmat(3,3),cmat(3*MAXICS)

  CHARACTER(LEN=2) :: cdir

  !*
  ! Start the exectuable code
  !!---------------------------

  !IF (dintv .EQ. 0.d0) THEN
  !  xsat(1:6)=SAT.x(1:6)
  !  RETURN
  !END IF

  CKF.mjd0=mjd
  CKF.sod0=sod
  CKF.mjd1 = INT(CKF.mjd0+(CKF.sod0+dintv)/86400.d0)
  CKF.sod1 = CKF.sod0+dintv-(CKF.mjd1-CKF.mjd0)*86400.d0
  CKF.rmjd=mjd
  CKF.rsod=sod

  ! decide if forward and/or backward integeration are needed
  sod0=(CKF.mjd0-CKF.rmjd)*86400.d0+CKF.sod0
  sod1=(CKF.mjd1-CKF.rmjd)*86400.d0+CKF.sod1

  tb=sod0
  te=sod1
       
  CKF.lgfm = .FALSE.
  CKF.ngc = 0
  CKF.ltog=0

  CKF.nprn=1
  CKF.cprn(1)=SAT.cprn

  CKF.lpart=.TRUE.

  CKF.nequ=6*(SAT.npar-6+1+6)

  CKF.int_step=15.d0

  nsign=SIGN(1.d0,dintv)

  IF (nsign .EQ. -1) THEN
    tdir=tb
  ELSE
    tdir=te
  END IF

  DO isat=1, CKF.nprn
    
    CALL ppp_fright_init(CKF,SAT)

    CALL ppp_fright_initorb(CKF.nequ,x,0)

    !! faster
    ! t=CKF.rsod
    ! hh=dintv
    ! !write(*,*) dintv
    ! CALL ppp_rkf(CKF.rmjd,t,x,hh,x,CKF.nequ,acc,amat,bmat,cmat)
    ! CALL ppp_fright_acc(CKF.rmjd,t,x,acc,amat,bmat,cmat)

    t=CKF.rsod
    hh=CKF.int_step*nsign
    ii=INT(ABS(dintv)/CKF.int_step)
    DO i=1, ii
     CALL ppp_rkf(CKF.rmjd,t,x,hh,x,CKF.nequ,acc,amat,bmat,cmat)
     t=t+hh
    END DO

    hh=dintv-ii*hh
    CALL ppp_rkf(CKF.rmjd,t,x,hh,x,CKF.nequ,acc,amat,bmat,cmat)

    ! t=t+hh
    ! CALL ppp_fright_acc(CKF.rmjd,t,x,acc,amat,bmat,cmat)

    ! Save state and transformation matrix
    ! km
    xsat(1:6)=x(1:6)
    SAT.phi(1:CKF.nequ-6)=x(7:CKF.nequ)
  END DO

  RETURN

END SUBROUTINE
