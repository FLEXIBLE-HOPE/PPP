!
SUBROUTINE dop(CKF,OB,SIT)
!!
!*
USE const
USE ckdctrl
USE station
USE observation
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!---------------------
TYPE(CKDCFG) :: CKF
TYPE(RNXOBS) :: OB
TYPE(SITE) :: SIT

  !*
  ! The local variables
  !!-----------------------
  INTEGER(IT) :: i,j,isat,ipar,isys
  INTEGER(IT) :: nsat(MAXSYS),lfn
  INTEGER(IT) :: iy,imon,id,ih,im
  REAL(RL) :: vdop(3),det,sec
  REAL(RL) :: norm(3,3),amat(3)

  LOGICAL lfirst
  DATA lfirst /.TRUE./
  SAVE lfirst,lfn

  !*
  ! The function called
  !!---------------------------
  INTEGER(IT) :: pointer_string
  INTEGER(IT) :: get_valid_unit

  !*
  ! Start the exectuable code
  !!---------------------------

  IF (pointer_string(OB.npar,OB.pname,'STAPX') .EQ. 0) RETURN

  IF (lfirst .EQ. .TRUE.) THEN
    lfirst=.FALSE.
    lfn=get_valid_unit(10)
    OPEN(UNIT=lfn,FILE='dop_'//TRIM(CKF.flnztd(5:)))
  END IF

  CALL mjd2date(CKF.mjd,CKF.sod,iy,imon,id,ih,im,sec)

  norm=0.d0
  nsat=0
  DO isat=1, CKF.nprn

    IF (OB.omc(isat,1) .EQ. 0.d0) CYCLE

    isys=INDEX(SYS,CKF.cprn(isat)(1:1))
    nsat(isys)=nsat(isys)+1

    amat=0.d0
    DO ipar=1,OB.npar
      SELECT CASE(TRIM(OB.pname(ipar)))
        CASE('STAPX')
          amat(1)=OB.amat(ipar,isat)
        CASE('STAPY')
          amat(2)=OB.amat(ipar,isat)
        CASE('STAPZ')
          amat(3)=OB.amat(ipar,isat)
      END SELECT
    END DO

    DO i=1, 3
      DO j=1, 3
        norm(i,j)=norm(i,j)+amat(i)*amat(j)
      END DO
    END DO

  END DO

  DO i=1, 3
    IF (norm(i,i) .EQ. 0.d0) norm(i,i)=norm(i,i)+1.d8
  END DO
  
  CALL matinv(norm,3,3,det)
  IF (det .EQ. 0.d0) RETURN

  vdop(1)=DSQRT(norm(1,1)+norm(2,2)+norm(3,3))
  vdop(2)=DSQRT(norm(1,1)+norm(2,2))
  vdop(3)=DSQRT(norm(3,3))
  SIT%dop(1:3)=vdop(1:3)

  ! PDOP, HDOP, VDOP, GPS, GLS, GAL, BDS, IRNSS, QZSS, SBAS, LEO
  WRITE(lfn,'(I4,4I3,F5.1,1X,A4,3F7.3,8I3)') iy,imon,id,ih,im,sec,SIT.name,vdop,nsat!,norm(1,1),norm(2,2)

  RETURN

END SUBROUTINE

!
SUBROUTINE dop_clk(CKF,OB,SIT)
!!
!*
USE const
USE ckdctrl
USE station
USE observation
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!---------------------
TYPE(CKDCFG) :: CKF
TYPE(RNXOBS) :: OB
TYPE(SITE) :: SIT

  !*
  ! The local variables
  !!-----------------------
  INTEGER(IT) :: i,j,isat,ipar,isys
  INTEGER(IT) :: nsat(MAXSYS),lfn
  INTEGER(IT) :: iy,imon,id,ih,im
  REAL(RL) :: vdop(3),det,sec
  !REAL(RL) :: norm(3+MAXSYS,3+MAXSYS),amat(3+MAXSYS)
  REAL(RL) :: norm(3+CKF.nsys,3+CKF.nsys),amat(3+MAXSYS)

  LOGICAL lfirst
  DATA lfirst /.TRUE./
  SAVE lfirst,lfn

  !*
  ! The function called
  !!---------------------------
  INTEGER(IT) :: pointer_string
  INTEGER(IT) :: get_valid_unit

  !*
  ! Start the exectuable code
  !!---------------------------

  IF (pointer_string(OB.npar,OB.pname,'STAPX') .EQ. 0) RETURN

  IF (lfirst .EQ. .TRUE.) THEN
    lfirst=.FALSE.
    lfn=get_valid_unit(10)
    OPEN(UNIT=lfn,FILE='dop_'//TRIM(CKF.flnztd(5:)))
  END IF

  CALL mjd2date(CKF.mjd,CKF.sod,iy,imon,id,ih,im,sec)

  norm=0.d0
  nsat=0
  DO isat=1, CKF.nprn

    IF (OB.omc(isat,1) .EQ. 0.d0) CYCLE

    isys=INDEX(SYS,CKF.cprn(isat)(1:1))
    nsat(isys)=nsat(isys)+1

    amat=0.d0
    DO ipar=1,OB.npar
      SELECT CASE(TRIM(OB.pname(ipar)))
        CASE('STAPX')
          amat(1)=OB.amat(ipar,isat)
        CASE('STAPY')
          amat(2)=OB.amat(ipar,isat)
        CASE('STAPZ')
          amat(3)=OB.amat(ipar,isat)
      END SELECT
    END DO

    isys=INDEX(CKF.system,CKF.cprn(isat)(1:1))
    amat(3+isys)=-1.d0
    DO i=1, 3+CKF.nsys
      DO j=1, 3+CKF.nsys
        norm(i,j)=norm(i,j)+amat(i)*amat(j)
      END DO
    END DO

  END DO

  DO i=1, 3+CKF.nsys
    IF (norm(i,i) .EQ. 0.d0) norm(i,i)=norm(i,i)+1.d8
  END DO
  
  CALL matinv(norm,3+CKF.nsys,3+CKF.nsys,det)
  IF (det .EQ. 0.d0) RETURN

  vdop(1)=DSQRT(norm(1,1)+norm(2,2)+norm(3,3))
  vdop(2)=DSQRT(norm(1,1)+norm(2,2))
  vdop(3)=DSQRT(norm(3,3))
  SIT.dop(1:3)=vdop(1:3)

  ! PDOP, HDOP, VDOP, GPS, GLS, GAL, BDS, IRNSS, QZSS, SBAS, LEO
  WRITE(lfn,'(I4,4I3,F5.1,1X,A4,3F7.3,8I3)') iy,imon,id,ih,im,sec,SIT.name,vdop,nsat

  RETURN

END SUBROUTINE
