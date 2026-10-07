!*
SUBROUTINE galileo_att(mjd,sod,xsat,xsun,xscf,yscf,zscf,cprn,blk)
!!
!*
USE const
USE ISO_FORTRAN_ENV
IMPLICIT NONE

  !*
  ! The local variables
  !!-------------------------

  INTEGER(IT) :: nepo,iepo
  INTEGER(IT) :: ierr,isat,lfn
  INTEGER(IT) :: iy,imon,id,ih,im,mjd0(MAXSAT),mjd,i
  CHARACTER(LEN=*) :: cprn,blk

  REAL(RL) :: nonyaw,inityaw(MAXSAT),modyaw,dt,per
  REAL(RL) :: det,xsat(6),xsun(6),beta,u,phi,isec,xu,sod0(MAXSAT),sod
  REAL(RL) :: xscf(1:*),yscf(1:*),zscf(1:*),yangle
  REAL(RL) SANTX, SANTY, v(3),r(3),s_unit(3)

  LOGICAL(LG) :: lnoon(MAXSAT),lmidnight(MAXSAT)
  DATA lnoon,lmidnight /.TRUE.,.TRUE./
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

  modyaw=DATAN2(-DTAN(beta),DSIN(u))*RAD2DEG
  nonyaw=DATAN2(-DTAN(beta),DSIN(u))*RAD2DEG

!!!! Locate the specific satellite according to the cprn.
  READ(cprn(2:3),'(I2)') isat

  xu=u
  IF (u .GT. PI) xu=u-2*PI

  IF (TRIM(blk) .EQ. 'GALILEO-2') THEN
    IF (DABS(beta)*RAD2DEG.LT.4.1d0 .AND. DABS(xu)*RAD2DEG.LT.10.d0) THEN
      IF (lmidnight(isat) .EQ. .TRUE.) THEN
        lmidnight(isat)=.FALSE.
        lnoon(isat)=.TRUE.
        inityaw(isat)=nonyaw
        mjd0(isat)=mjd
        sod0(isat)=sod
!        beta0=beta
        WRITE(*,'(I5,4I3,F11.7,1X,A3,1X,A20,5F8.2)') iy,imon,id,ih,im,isec,cprn,TRIM(blk)//'*',beta*RAD2DEG,u*RAD2DEG,inityaw(isat),SIGN(1.d0,inityaw(isat))
      END IF
      dt=timdif(mjd,sod,mjd0(isat),sod0(isat))
      modyaw=(90*SIGN(1.d0,inityaw(isat))+(inityaw(isat)-90*SIGN(1.d0,inityaw(isat)))*DCOS(2*PI/5656.d0*dt))
    END IF
  
    IF (DABS(beta)*RAD2DEG.LT.4.1d0 .AND. DABS(u-PI)*RAD2DEG.LT.10.d0) THEN
      IF (lnoon(isat) .EQ. .TRUE.) THEN
        lmidnight(isat)=.TRUE.
        lnoon(isat)=.FALSE.
        inityaw(isat)=nonyaw
        mjd0(isat)=mjd
        sod0(isat)=sod
        WRITE(*,'(I5,4I3,F11.7,1X,A3,1X,A20,5F8.2)') iy,imon,id,ih,im,isec,cprn,TRIM(blk),beta*RAD2DEG,u*RAD2DEG,inityaw(isat),SIGN(1.d0,inityaw(isat))
      END IF
      dt=timdif(mjd,sod,mjd0(isat),sod0(isat))
      modyaw=(90*SIGN(1.d0,inityaw(isat))+(inityaw(isat)-90*SIGN(1.d0,inityaw(isat)))*DCOS(2*PI/5656.d0*dt))
    END IF
  ELSE IF (TRIM(blk) .EQ. 'GALILEO-1') THEN
        s_unit(1)=-DSIN(u-PI)*DCOS(beta)
        s_unit(2)=-DSIN(beta)
        s_unit(3)=-DCOS(u-PI)*DCOS(beta)
        IF (DABS(s_unit(1)).LT.DSIN(15.d0*DEG2RAD) .AND. DABS(s_unit(2)).LT.DSIN(2.d0*DEG2RAD)) THEN
          xu=u
          IF (u .GT. PI) xu=u-2*PI
          IF (lmidnight(isat).EQ..TRUE. .AND. DABS(xu)*RAD2DEG.LE.30) THEN
            lmidnight(isat)=.FALSE.
            lnoon(isat)=.TRUE.
            inityaw(isat)=s_unit(2)
            WRITE(*,'(I5,4I3,F11.7,1X,A3,1X,A20,5F8.2)') iy,im,id,ih,im,isec,cprn,TRIM(blk)//'*',beta*RAD2DEG,u*RAD2DEG,inityaw(isat),SIGN(1.d0,inityaw(isat))
          END IF

          IF (lnoon(isat).EQ..TRUE. .AND. DABS(u-PI)*RAD2DEG.LE.30) THEN
            lnoon(isat)=.FALSE.
            lmidnight(isat)=.TRUE.
            inityaw(isat)=s_unit(2)
            WRITE(*,'(I5,4I3,F11.7,1X,A3,1X,A20,5F8.2)') iy,im,id,ih,im,isec,cprn,TRIM(blk),beta*RAD2DEG,u*RAD2DEG,inityaw(isat),SIGN(1.d0,inityaw(isat))
          END IF
          s_unit(2)=0.5d0*(DSIN(2.d0*DEG2RAD)*SIGN(1.d0,inityaw(isat))+s_unit(2))+0.5d0*(DSIN(2.d0*DEG2RAD)*SIGN(1.d0,inityaw(isat))-s_unit(2))*DCOS(PI*DABS(s_unit(1))/DSIN(15.d0*DEG2RAD))
          s_unit(3)=DSQRT(1-s_unit(1)**2-s_unit(2)**2)*SIGN(1.d0,s_unit(3))
          nonyaw=DATAN2(-s_unit(2)/DSQRT(1-s_unit(3)**2),-s_unit(1)/DSQRT(1-s_unit(3)**2))*RAD2DEG
        END IF
  END IF

  DO i=1,3
       R(i)=xsat(i)/DSQRT(xsat(1)**2+xsat(2)**2+xsat(3)**2)
       V(i)=xsat(i+3)/DSQRT(xsat(4)**2+xsat(5)**2+xsat(6)**2)
  END DO

  YANGLE=DACOS((xscf(1)*xsat(4)+xscf(2)*xsat(5)+xscf(3)*xsat(6))/ &
                DSQRT(xsat(4)**2+xsat(5)**2+xsat(6)**2))*RAD2DEG
  IF (beta .GT. 0.d0) yangle=-yangle

  SANTX=(DCOS((modyaw-YANGLE)*DEG2RAD)*(V(2)-V(3)*R(2)/R(3))-DCOS(modyaw*DEG2RAD)* &
         (xscf(2)-xscf(3)*R(2)/R(3)))/(xscf(1)*V(2)-xscf(2)*v(1)+ &
         ((xscf(2)*V(3)-xscf(3)*V(2))*R(1)+(xscf(3)*V(1)- &
         xscf(1)*V(3))*R(2))/R(3))
  SANTY=(DCOS(modyaw*DEG2RAD) - (V(1)-V(3)*R(1)/R(3))*SANTX)/(V(2)-V(3)*R(2)/R(3))

  ! THE BODY-X UNIT VECTOR ROTATED BY (modyaw-YANGLE) RETURNED
  IF (TRIM(blk) .EQ. 'GALILEO-2') THEN
    IF ( (DABS(beta)*RAD2DEG.LT.4.1d0 .AND. DABS(xu)*RAD2DEG.LT.10.d0) .OR. (DABS(beta)*RAD2DEG.LT.4.1d0 .AND. DABS(u-PI)*RAD2DEG.LT.10.d0) ) THEN
      xscf(1)=SANTX
      xscf(2)=SANTY
      xscf(3)=(-R(1)*SANTX-R(2)*SANTY)/R(3)

      CALL CROSS(zscf,xscf,yscf)
      CALL UNIT_VECTOR(3,yscf,yscf,DET)
!    WRITE(*,*) mjd,sod,mjd0(isat),sod0(isat),modyaw,yangle
    END IF 
  ELSE IF(TRIM(blk) .EQ. 'GALILEO-1' ) THEN
    IF (DABS(s_unit(1)).LT.DSIN(15.d0*DEG2RAD) .AND. DABS(s_unit(2)).LT.DSIN(2.d0*DEG2RAD)) THEN
      xscf(1)=SANTX
      xscf(2)=SANTY
      xscf(3)=(-R(1)*SANTX-R(2)*SANTY)/R(3)

      CALL CROSS(zscf,xscf,yscf)
      CALL UNIT_VECTOR(3,yscf,yscf,DET)
!    WRITE(*,*) mjd,sod,mjd0(isat),sod0(isat),modyaw,yangle
    END IF 
  END IF


END SUBROUTINE
