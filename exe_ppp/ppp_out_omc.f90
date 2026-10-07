!*
!!输出omc
!!xsy:2022/8/29
!*
SUBROUTINE ppp_out_omc(CKF,SIT,OB,SAT)
!!
!*
USE ckdctrl
USE observation
USE station
USE satellite 
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!*
! The arguments
!!---------------------
TYPE(CKDCFG) :: CKF
TYPE(RNXOBS) :: OB
TYPE(SITE) :: SIT
TYPE(SATE) :: SAT(1:*)
   
  !*
  ! The local variables
  !!-----------------------
  INTEGER(IT) :: i,isit,isat,iy,imon,id,ih,im,outsat,numnlos,ipar,indstapx
  REAL(RL) :: sec,omc_if
  
   outsat = 0
   numnlos = 0
  DO isat=1,CKF.nprn
    IF (OB.omc(isat,1+MAXFREQ) .NE. 0.d0) THEN
      outsat = outsat + 1
      IF (OB.nlosflag(isat) .EQ. 0) THEN
        numnlos = numnlos + 1
      END IF
    END IF
  END DO

  CALL mjd2date(CKF.mjd,CKF.sod,iy,imon,id,ih,im,sec)
  WRITE(1003,'(A3,I5,4I3,F11.7,I7,F10.2,1X,I3,1X,I3,4X,A)') 'TIM',iy,imon,id,ih,im,sec,CKF.mjd,CKF.sod,numnlos,outsat,'Valid_Range_OMC'
  ! WRITE(5000,'(A)') '#---MJD-----SOD--SIT-PRN-----------X0-----------Y0-----------Z0------RCLK0-------ZTD0-------ION0-------OMC_P1------OMC_P2------OMC_L1------OMC_L2-----ZTDMAP-------XPAR-------YPAR-------ZPAR---CSLP---ELEV'
  ! ipar=0
  ! DO WHILE(ipar .LT. OB%npar)
  !   ipar=ipar+1
  !   IF (OB%pname(ipar)(1:5) .EQ. 'STAPX') THEN
  !     indstapx=ipar
  !     EXIT
  !   END IF
  ! END DO
  DO isat=1, CKF.nprn
    IF (OB.omc(isat,MAXFREQ+1).NE.0.d0 .AND. OB.omc(isat,MAXFREQ+2).NE.0.d0) THEN
      omc_if = OB.omc(isat,MAXFREQ+1)*SAT(isat).fac(1)-OB.omc(isat,MAXFREQ+2)*SAT(isat).fac(2)
      WRITE(1003,'(2X,A4,A4,I3,5F16.8)')SIT.name,CKF.cprn(isat),OB.nlosflag(isat),OB.omc(isat,MAXFREQ+1),OB.omc(isat,MAXFREQ+2),OB.omc(isat,MAXFREQ+3),omc_if,OB.elev(isat)*RAD2DEG
    ELSE
      omc_if = 0.d0
      IF (OB.elev(isat).GE.SIT%cutoff) THEN
        WRITE(1003,'(2X,A4,A4,I3,5F16.8)')SIT.name,CKF.cprn(isat),OB.nlosflag(isat),OB.omc(isat,MAXFREQ+1),OB.omc(isat,MAXFREQ+2),OB.omc(isat,MAXFREQ+3),omc_if,OB.elev(isat)*RAD2DEG
      END IF
    END IF

    ! IF (OB%omc(isat,MAXFREQ+1).NE.0.D0) THEN
    !   WRITE(5000,'(I7,F8.1,1X,A4,1X,A3,1X,3(F12.3,1X),3(F10.3,1X),4(F12.3),4(1X,F10.3),1X,I6,1X,F6.1)') CKF%mjd,CKF%sod,SIT%name,CKF%cprn(isat),&
    !                 SIT%x(1:3),SIT%rclock(CKF%iref),&
    !                 SIT%trpdel(isat),SIT%iondel(isat,1),&
    !                 OB%omc(isat,MAXFREQ+1),OB%omc(isat,MAXFREQ+2),OB%omc(isat,1),OB%omc(isat,2),&
    !                 OB%zmap(isat),OB%amat(indstapx:indstapx+2,isat),OB%flag(isat,1),OB%elev(isat)*RAD2DEG
    ! END IF
  END DO
  ! WRITE(5000,'(A)') '#-PNAME------------XINI---------XCOR---------XEST----MAP-----------RW'

  RETURN
  
END SUBROUTINE
