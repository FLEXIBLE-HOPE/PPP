!*
SUBROUTINE ppp_spp(CKF,SAT,SIT,OB,corx)
!!
!*
USE ckdctrl
USE station
USE satellite
USE observation
IMPLICIT NONE

!*
! The arguments
!!--------------------------
TYPE(CKDCFG) :: CKF
TYPE(SATE) :: SAT(MAXSAT)
TYPE(SITE) :: SIT
TYPE(RNXOBS) :: OB
REAL(RL) :: corx(5)

  !*
  ! The local variables
  !!------------------------------
  INTEGER(IT) :: i,j,k,ipar
  INTEGER(IT) :: isat,isys,nobs,nos(MAXSYS)

  REAL(RL) :: omc(MAXSAT),clk(MAXSYS)
  REAL(RL) :: amat(MAXSYS+3,MAXSAT)
  REAL(RL) :: norm(4,4)
  REAL(RL) :: wl(MAXSYS+3),delt(MAXSYS+3)

  !*
  ! The function called
  !!------------------------------
  INTEGER(IT) :: pointer_string

  !*
  ! Start the exectuable code
  !!------------------------------

  corx=0.d0

  ipar=0
  ipar=pointer_string(OB.npar,OB.pname,'STAPX')
  IF (ipar .EQ. 0) THEN
    WRITE(OUTPUT_UNIT,'(A)') '***WARNING(lsq_spp): kinematic station without coordinates to be estimated'
    RETURN
  END IF

  nos=0
  clk=0.d0
  DO isat=1, CKF.nprn
    isys=INDEX(SYS,CKF.cprn(isat)(1:1))
    IF (OB.omc(isat,1).NE.0.d0 .AND. OB.omc(isat,MAXFREQ+1).NE.0.d0) THEN
      nos(isys)=nos(isys)+1
      clk(isys)=clk(isys)+OB.omc(isat,MAXFREQ+1)*SAT(isat).fac(1)-OB.omc(isat,MAXFREQ+2)*SAT(isat).fac(2)
    END IF
  END DO
  DO isys=1, MAXSYS
    IF (nos(isys) .NE. 0) THEN
      clk(isys)=clk(isys)/nos(isys)
    END IF
  END DO


  nobs=0
  amat=0.d0
  omc=0.d0
  DO isat=1, CKF.nprn
    isys=INDEX(SYS,CKF.cprn(isat)(1:1))

    IF (OB.omc(isat,1).NE.0.d0 .AND. OB.omc(isat,MAXFREQ+1).NE.0.d0) THEN
      nobs=nobs+1
      omc(nobs)=OB.omc(isat,MAXFREQ+1)*SAT(isat).fac(1)-OB.omc(isat,MAXFREQ+2)*SAT(isat).fac(2)-clk(isys)

      ipar=0
      DO WHILE(ipar .LT. OB.npar)
        ipar=ipar+1
        IF (OB.ltog(ipar,isat) .EQ. 0) CYCLE
        IF (OB.pname(ipar)(1:5) .EQ. 'STAPX') THEN
          DO i=1,3
            amat(i,nobs)=OB.amat(ipar+i-1,isat)
          END DO
          ipar=ipar+2
        ELSE IF (OB.pname(ipar)(1:6) .EQ. 'RECCLK') THEN
          amat(4,nobs)=1.d0
        END IF
      END DO
    END IF
  END DO

  norm=0.d0
  wl=0.d0
  delt=0.d0
  IF (nobs .GT. 4) THEN
    DO i=1, 4
      DO j=1, nobs
        wl(i)=wl(i)+amat(i,j)*omc(j)
      END DO
      DO j=1, i
        DO k=1, nobs
          norm(i,j)=norm(i,j)+amat(i,k)*amat(j,k)
        END DO
        norm(j,i)=norm(i,j)
      END DO
    END DO

    DO i=1, 4
      IF (norm(i,i) .EQ. 0.d0) norm(i,i)=norm(i,i)+1.d8
    END DO

    CALL matinv(norm,4,4,delt(1))
    IF (delt(1) .EQ. 0.d0) RETURN

    DO i=1, 4
      DO j=1, 4
        delt(i)=delt(i)+norm(i,j)*wl(j)
      END DO
    END DO

    corx(5)=0.d0
    DO i=1, 3
      corx(i)=delt(i)
      corx(5)=corx(5)+delt(i)*delt(i)
    END DO
    corx(4)=delt(4)
    corx(5)=DSQRT(corx(5))
  ELSE
    RETURN
  END IF

  RETURN

END SUBROUTINE
