!*
SUBROUTINE beidou3_att(mjd,sod,xsat,xsun,xscf,yscf,zscf,cprn,blk)
!!
!! Bug need to fix if the first epoch is already in the requested period, which could result in a wrong indication of the yaw start epoch.
!! Similar epoch file like bds_attribute_epoch needed?
!*
USE const
USE ISO_FORTRAN_ENV
IMPLICIT NONE

  !*
  ! The local variables
  !!-------------------------

  INTEGER(IT) :: nepo,iepo
  INTEGER(IT) :: ierr,isat,lfn
  INTEGER(IT) :: iy,imon,id,ih,im,mjd0(MAXSAT),mjd,i,mjds
  CHARACTER(LEN=*) :: cprn,blk

  REAL(RL) :: nonyaw,inityaw(MAXSAT),modyaw,dt,per
  REAL(RL) :: det,xsat(6),xsun(6),beta,u,phi,isec,xu,sod0(MAXSAT),sod,sods
  REAL(RL) :: xscf(1:*),yscf(1:*),zscf(1:*),yangle
  REAL(RL) SANTX, SANTY, v(3),r(3)
  REAL(RL) :: f,betad,phid

  LOGICAL(LG) :: lnoon(MAXSAT),lmidnight(MAXSAT),lfirst
  DATA lnoon,lmidnight /.TRUE.,.TRUE./
  DATA lfirst /.TRUE./
  SAVE mjd0,sod0,lnoon,lmidnight,inityaw

  !*
  ! The function called
  !!--------------------------
  INTEGER(IT) :: get_valid_unit
  REAL(RL) :: dot,timdif

  !*
  ! Start the exectuable code
  !!--------------------------


  CALL betau(xsat,xsun,beta,u)

  xu=u
  IF (u .GT. PI) xu=u-2*PI
  CALL mjd2date(mjd,sod,iy,imon,id,ih,im,isec)

  modyaw=DATAN2(-DTAN(beta),DSIN(u))*RAD2DEG
  nonyaw=DATAN2(-DTAN(beta),DSIN(u))*RAD2DEG

  f=1.d0/(1.d0+80000.d0*(DSIN(u)**4))
  IF (blk(1:10) .EQ. 'BEIDOU-3MC'.OR. blk(1:10) .EQ. 'BEIDOU-3IC') THEN
     betad=beta+f*(SIGN(3.d0*DEG2RAD,beta)-beta)
  ELSE IF (blk(1:10) .EQ. 'BEIDOU-3MS' .OR. blk(1:10) .EQ. 'BEIDOU-3IS' ) THEN
     betad=3.d0*DEG2RAD*SIGN(1.d0,beta)
  ELSE
     betad=beta+f*(SIGN(3.d0*DEG2RAD,beta)-beta)
  END IF
  !WRITE(*,'(10F8.2)') f,betad,beta

  phid=0.d0
  IF (DABS(beta)*RAD2DEG.LT. 3.d0 ) THEN
     phid=DATAN2(-DTAN(betad),DSIN(u))
  ELSE
     phid=DATAN2(-DTAN(beta),DSIN(u))
  END IF
  phid=phid*RAD2DEG

!  WRITE(3000,'(I5,4I3,F11.7,1X,A3,1X,A20,10F8.2)') iy,imon,id,ih,im,isec,cprn,blk,beta*RAD2DEG,u*RAD2DEG,nonyaw,nonyaw,0.d0,phid

  DO i=1,3
       R(i)=xsat(i)/DSQRT(xsat(1)**2+xsat(2)**2+xsat(3)**2)
       V(i)=xsat(i+3)/DSQRT(xsat(4)**2+xsat(5)**2+xsat(6)**2)
  END DO

  YANGLE=DACOS((xscf(1)*xsat(4)+xscf(2)*xsat(5)+xscf(3)*xsat(6))/ &
                DSQRT(xsat(4)**2+xsat(5)**2+xsat(6)**2))*RAD2DEG
  IF (beta .GT. 0.d0) yangle=-yangle

  SANTX=(DCOS((phid-YANGLE)*DEG2RAD)*(V(2)-V(3)*R(2)/R(3))-DCOS(phid*DEG2RAD)* &
         (xscf(2)-xscf(3)*R(2)/R(3)))/(xscf(1)*V(2)-xscf(2)*v(1)+ &
         ((xscf(2)*V(3)-xscf(3)*V(2))*R(1)+(xscf(3)*V(1)- &
         xscf(1)*V(3))*R(2))/R(3))
  SANTY=(DCOS(phid*DEG2RAD) - (V(1)-V(3)*R(1)/R(3))*SANTX)/(V(2)-V(3)*R(2)/R(3))

  ! THE BODY-X UNIT VECTOR ROTATED BY (modyaw-YANGLE) RETURNED
  IF ( (DABS(beta)*RAD2DEG.LT.3.d0 .AND. DABS(xu)*RAD2DEG.LT.10.d0) .OR. (DABS(beta)*RAD2DEG.LT.3.d0 .AND. DABS(u-PI)*RAD2DEG.LT.10.d0) ) THEN
    xscf(1)=SANTX
    xscf(2)=SANTY
    xscf(3)=(-R(1)*SANTX-R(2)*SANTY)/R(3)

    CALL CROSS(zscf,xscf,yscf)
    CALL UNIT_VECTOR(3,yscf,yscf,DET)
!    WRITE(*,*) mjd,sod,mjd0(isat),sod0(isat),modyaw,yangle
  END IF


END SUBROUTINE
