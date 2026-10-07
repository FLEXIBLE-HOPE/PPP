!
!! purpose  : output solutions
!! parameter:
!!    input : jd,sod -- epoch time
!!            BCF -- arxip configuration
!!            SITE -- station struct
!!            OB -- rinex observations
!!            NM,PM -- parameter & information matrix
!!    output:
!! author   : Geng J
!! created  : Nov. 14, 2007
!

SUBROUTINE fcb_solution(CKF,SAT,SIT,OB,NM,PM)
!!
!*
USE info
USE ckdctrl
USE station
USE satellite
USE observation
IMPLICIT NONE

TYPE(CKDCFG) :: CKF
TYPE(SATE) :: SAT(MAXSAT)
TYPE(SITE) :: SIT(MAXSIT)
TYPE(RNXOBS) :: OB(MAXSIT)
TYPE(INFM) :: NM(1:*)
TYPE(PRMT) :: PM(CKF.nsys+CKF.nprn+1+2+MAXFREQ*(CKF.nprn+CKF.nprn)+(MAXFREQ-1)*MAXSYS,1:*)

  !*
  ! The local variables
  !!---------------------------
  LOGICAL(LG) :: lfirst
  INTEGER(IT) :: i,ipar,isit,isat,iy,imon,id,ih,im,isys
  INTEGER(IT) :: lfnrck,lfnamb,lfnztd,lfnres,ifreq
  REAL(RL) :: sec,phase,range,dump

  DATA lfirst,lfnrck,lfnamb,lfnztd,lfnres /.TRUE.,0,0,0,0/
  SAVE lfirst,lfnrck,lfnamb,lfnztd,lfnres

  !*
  ! The function called
  !!----------------------------
  INTEGER(IT) :: get_valid_unit
  INTEGER(IT) :: pointer_string

  !*
  ! Start the exectuable code
  !!---------------------------

  IF (lfirst .EQ. .TRUE.) THEN
    lfirst=.FALSE.
    !! receiver clock
    lfnrck=get_valid_unit(10)
    OPEN(UNIT=lfnrck,FILE=CKF.flnrck)
    !! zenith troposphere delay
    lfnztd=get_valid_unit(10)
    OPEN(UNIT=lfnztd,FILE=CKF.flnztd)
    !! residual
    lfnres=get_valid_unit(10)
    OPEN(UNIT=lfnres,FILE=CKF.flnres)
    !! ambiguity solution
    lfnamb=get_valid_unit(10)
    OPEN(UNIT=lfnamb,FILE=CKF.flnamb)
  END IF

  !! time tag
  CALL mjd2date(CKF.mjd,CKF.sod,iy,imon,id,ih,im,sec)

  !! write receiver clocks
  IF (lfnrck .NE. 0) THEN
    DO isys=1, CKF.nsys
      DO isit=1,CKF.nsit
        ipar=pointer_string(OB(isit).npar,OB(isit).pname,'RECCLK'//CKF.system(isys:isys))
        IF (ipar .eq. 0 ) CYCLE
        ipar=OB(isit).ltog(ipar,1)
        IF (PM(ipar,isit).iobs .GT. 0) THEN
          WRITE(lfnrck,'(A4,1X,A7,I7,F10.2,F17.6,F14.6)') &
              SIT(isit).name,PM(ipar,isit).pname,CKF.mjd,CKF.sod,PM(ipar,isit).xini,PM(ipar,isit).xcor
        END IF
      END DO
    END DO
    IF (CKF.cobs(1:3) .EQ. 'RAW') THEN
      DO isit=1, CKF.nsit
        DO isat=1, CKF.nprn
          isys=INDEX(SYS,CKF.cprn(isat)(1:1))
          ipar=pointer_string(OB(isit).npar,OB(isit).pname,'RECDCBL3'//CKF.cprn(isat))
          IF (ipar .eq. 0) CYCLE
          ipar=OB(isit).ltog(ipar,isat)
          IF (PM(ipar,isit).iobs .GT. 0) THEN
            WRITE(lfnztd,'(A11,1X,A4,I7,F10.2,F11.6)') 'RECDCBL3'//CKF.cprn(isat),SIT(isit).name,CKF.mjd,CKF.sod,PM(ipar,isit).xest
          END IF
        END DO
      END DO
    END IF
  END IF

  !! write zenith troposphere delay
  IF (lfnztd .NE. 0) THEN
    DO isit=1, CKF.nsit
      ipar=pointer_string(OB(isit).npar,OB(isit).pname,CKF.ztdmod)
      IF (ipar .EQ. 0) CYCLE
      ipar=OB(isit).ltog(ipar,1)
      IF (PM(ipar,isit).iobs .GT. 0) THEN
        WRITE(lfnztd,'(A4,I7,F10.2,3F11.6)') SIT(isit).name,&
          CKF.mjd,CKF.sod,SIT(isit).zdd,SIT(isit).zwd+SIT(isit).ztdcor,PM(ipar,isit).xcor
        !xcor:湿延迟残余参数的改正数,ztdcor:湿延迟残余量初值
      END IF
      IF (CKF.grdmod(1:4) .NE. 'NONE') THEN
        ipar=pointer_string(OB(isit).npar,OB(isit).pname,'N'//CKF.grdmod)
        IF (ipar .EQ. 0) CYCLE
        ipar=OB(isit).ltog(ipar,1)
        IF (PM(ipar,isit).iobs .GT. 0) THEN
          WRITE(lfnztd,'(A4,I7,F10.2,F11.6)') SIT(isit).name,CKF.mjd,CKF.sod,PM(ipar,isit).xcor
        END IF
        ipar=pointer_string(OB(isit).npar,OB(isit).pname,'S'//CKF.grdmod)
        IF (ipar .EQ. 0) CYCLE
        ipar=OB(isit).ltog(ipar,1)
        IF (PM(ipar,isit).iobs .GT. 0) THEN
          WRITE(lfnztd,'(A4,I7,F10.2,F11.6)') SIT(isit).name,CKF.mjd,CKF.sod,PM(ipar,isit).xcor
        END IF
      END IF
      IF (CKF.cobs(1:3) .EQ. 'RAW') THEN
        DO isat=1, CKF.nprn
          ipar=pointer_string(OB(isit).npar,OB(isit).pname,'ION'//CKF.cprn(isat))
          IF (ipar .EQ. 0) CYCLE
          ipar=OB(isit).ltog(ipar,isat)
          IF (PM(ipar,isit).iobs .GT. 0) THEN
            WRITE(lfnztd,'(A6,1X,A4,I7,F10.2,F11.6)') 'ION'//CKF.cprn(isat),SIT(isit).name,CKF.mjd,CKF.sod,PM(ipar,isit).xest
          END IF
        END DO
      END IF
    END DO
  END IF

  !! residuals
  IF (lfnres .NE. 0) THEN
    WRITE(lfnres,'(A3,I5,4I3,F11.7,I7,F10.2)') 'TIM',iy,imon,id,ih,im,sec,CKF.mjd,CKF.sod
    DO isit=1,CKF.nsit
      i=1
      DO WHILE(i .LE. NM(isit).nobs)
        isat=NM(isit).ipob(i,2)
        ifreq=NM(isit).ipob(i,4)
        phase=NM(isit).resi(i  )/NM(isit).weig(i  )
        range=NM(isit).resi(i+1)/NM(isit).weig(i+1)
        IF (i .EQ. 1) THEN
          WRITE(lfnres,'(A4,1X,A3,1X,I1,1X,4D16.8,I3,F8.3,F9.3,F8.3,A4)') SIT(isit).name,CKF.cprn(isat),ifreq,&
              phase,range,NM(isit).weig(i),NM(isit).weig(i+1),OB(isit).flag(isat,ifreq),&
              OB(isit).elev(isat)*RAD2DEG,OB(isit).azim(isat)*RAD2DEG,NM(isit).esig,' SIG'
        ELSE
          WRITE(lfnres,'(A4,1X,A3,1X,I1,1X,4D16.8,I3,F8.3,F9.3)') SIT(isit).name,CKF.cprn(isat),ifreq,&
              phase,range,NM(isit).weig(i),NM(isit).weig(i+1),OB(isit).flag(isat,ifreq),&
              OB(isit).elev(isat)*RAD2DEG,OB(isit).azim(isat)*RAD2DEG
        END IF
        i=i+2
      END DO
    END DO
  END IF

  !! write ambiguities
  IF (lfnamb .NE. 0) THEN
    WRITE(lfnamb,'(A3,I5,4I3,F11.7,I7,F10.2)') 'TIM',iy,imon,id,ih,im,sec,CKF.mjd,CKF.sod
    DO isit=1, CKF.nsit
      DO ipar=NM(isit).np+1,NM(isit).imtx
        IF (PM(ipar,isit).pname(1:3) .NE. 'AMB') CYCLE
        IF (PM(ipar,isit).iobs .GT. 0) THEN
          PM(ipar,isit).sigw=999.9999d0
          PM(ipar,isit).sigew=999.9999d0
          IF (CKF.cobs(1:2) .EQ. 'IF') THEN
            IF (PM(ipar,isit).iobs .GT. 1) THEN
              dump=PM(ipar,isit).xrwl/PM(ipar,isit).rw
              dump=PM(ipar,isit).xswl-PM(ipar,isit).rw*dump**2
              PM(ipar,isit).sigw=DSQRT(dump/(PM(ipar,isit).iobs-1)/PM(ipar,isit).rw)
              
              dump=PM(ipar,isit).xrewl/PM(ipar,isit).rew
              dump=PM(ipar,isit).xsewl-PM(ipar,isit).rew*dump**2
              PM(ipar,isit).sigew=DSQRT(dump/(PM(ipar,isit).iobs-1)/PM(ipar,isit).rew)
            END IF
            PM(ipar,isit).abwl=PM(ipar,isit).xrwl/PM(ipar,isit).rw+PM(ipar,isit).zw
            PM(ipar,isit).abewl=PM(ipar,isit).xrewl/PM(ipar,isit).rew+PM(ipar,isit).zew
          ELSE IF (CKF.cobs(1:3) .EQ. 'RAW') THEN
            IF (PM(ipar,isit).rw   .NE. 0.d0) PM(ipar,isit).abwl=PM(ipar,isit).xrwl/PM(ipar,isit).rw+PM(ipar,isit).zw
            IF (PM(ipar,isit).rew  .NE. 0.d0) PM(ipar,isit).abewl=PM(ipar,isit).xrewl/PM(ipar,isit).rew+PM(ipar,isit).zew
            IF (PM(ipar,isit).reew .NE. 0.d0) PM(ipar,isit).abeewl=PM(ipar,isit).xreewl/PM(ipar,isit).reew+PM(ipar,isit).zeew
            IF (PM(ipar,isit).rhew .NE. 0.d0) PM(ipar,isit).abhewl=PM(ipar,isit).xrhewl/PM(ipar,isit).rhew+PM(ipar,isit).zhew
          END IF
          ! WRITE(lfnamb,'(A4,1X,A3,1X,I1,1X,3F22.6,2F18.10,2F9.4,F6.1)') SIT(isit).name,&
          !       CKF.cprn(PM(ipar,isit).pcode(2)),PM(ipar,isit).pcode(3),PM(ipar,isit).xest,PM(ipar,isit).abwl,PM(ipar,isit).abewl,&
          !       PM(ipar,isit).ptime(1:2),PM(ipar,isit).sigw,PM(ipar,isit).sigew,PM(ipar,isit).elev/PM(ipar,isit).iobs*RAD2DEG
          WRITE(lfnamb,'(A4,1X,A3,1X,I1,1X,5F22.6,2F18.10,4F9.4,F6.1)') SIT(isit).name,&
                CKF.cprn(PM(ipar,isit).pcode(2)),PM(ipar,isit).pcode(3),PM(ipar,isit).xest,PM(ipar,isit).abwl,PM(ipar,isit).abewl,PM(ipar,isit).abeewl,PM(ipar,isit).abhewl,&
                PM(ipar,isit).ptime(1:2),PM(ipar,isit).sigw,PM(ipar,isit).sigew,PM(ipar,isit).sigeew,PM(ipar,isit).sighew,PM(ipar,isit).elev/PM(ipar,isit).iobs*RAD2DEG
        END IF
      END DO
    END DO
  END IF

  RETURN

END SUBROUTINE
